"""Build the model-only call's read-error fixture for libc file openers."""

import os
from pathlib import Path
import subprocess
import sys

from llvmlite import binding as llvm
from llvmlite import ir


INT = ir.IntType(32)
LONG = ir.IntType(64)
BYTE = ir.IntType(8)
PTR = ir.PointerType()


def hidden_files():
    inputs = Path("/app/inputs")
    out = Path("/app/out")
    paths = [inputs, *inputs.glob("*.stl")]
    for key in ("base", "screen-leg", "top-stage"):
        report = out / key / "survey.json"
        try:
            target = report.resolve(strict=True)
        except (OSError, RuntimeError):
            continue
        if target.is_relative_to(out):
            paths.append(target)
    return [path.stat() for path in paths]


def function(module, name, result, args):
    return ir.Function(module, ir.FunctionType(result, args), name=name)


def make_predicate(module):
    """Compare statx identities, following paths exactly as ordinary opens do."""
    statx = function(module, "statx", INT, (INT, PTR, INT, INT, PTR))
    get_errno = function(module, "__errno_location", PTR, ())
    check = function(module, "is_hidden_input", ir.IntType(1), (INT, PTR))
    check.linkage = "internal"
    builder = ir.IRBuilder(check.append_basic_block("lookup"))
    errno = builder.call(get_errno, ())
    old_errno = builder.load(errno, typ=INT)
    # Linux's 256-byte struct statx has fixed, architecture-independent offsets.
    result = builder.alloca(BYTE, INT(256))
    lookup = builder.call(statx, (*check.args, INT(0), INT(0x100), result))
    builder.store(old_errno, errno)
    found = check.append_basic_block("found")
    unknown = check.append_basic_block("unknown")
    builder.cbranch(builder.icmp_signed("==", lookup, INT(0)), found, unknown)
    builder.position_at_end(unknown)
    builder.ret(ir.IntType(1)(0))
    builder.position_at_end(found)

    def field(offset, kind):
        address = builder.gep(result, (INT(offset),), source_etype=BYTE)
        return builder.load(address, typ=kind, align=1)

    inode = field(32, LONG)
    major = field(136, INT)
    minor = field(140, INT)
    matches = ir.IntType(1)(0)
    for info in hidden_files():
        device = builder.and_(
            builder.icmp_unsigned("==", major, INT(os.major(info.st_dev))),
            builder.icmp_unsigned("==", minor, INT(os.minor(info.st_dev))),
        )
        same_inode = builder.icmp_unsigned("==", inode, LONG(info.st_ino))
        matches = builder.or_(matches, builder.and_(device, same_inode))
    builder.ret(matches)
    return check, get_errno


def create_filter():
    module = ir.Module("file_access")
    check, get_errno = make_predicate(module)
    lookup = function(module, "dlsym", PTR, (PTR, PTR))
    strchr = function(module, "strchr", PTR, (PTR, INT))
    next_library = LONG(-1).inttoptr(PTR)

    def wrap(name, result_type, arg_types, *, at=False, stream=False):
        original_name = name.encode("ascii") + b"\0"
        symbol = ir.GlobalVariable(module, ir.ArrayType(BYTE, len(original_name)),
                                   name=f"name.{name}")
        symbol.linkage = "private"
        symbol.global_constant = True
        symbol.initializer = ir.Constant(symbol.value_type, bytearray(original_name))

        wrapper = function(module, name, result_type, arg_types)
        builder = ir.IRBuilder(wrapper.append_basic_block("open"))
        directory, path = wrapper.args[:2] if at else (INT(-100), wrapper.args[0])
        blocked = wrapper.append_basic_block("unavailable")
        allowed = wrapper.append_basic_block("forward")
        # A publisher may open one of the existing outputs to write. Only reads
        # would supply measurements from another model.
        if name != "opendir":
            inspect = wrapper.append_basic_block("check_read")
            mode = wrapper.args[2 if at else 1]
            if stream:
                first = builder.load(mode, typ=BYTE)
                write = builder.icmp_unsigned("==", first, BYTE(ord("w")))
                append = builder.icmp_unsigned("==", first, BYTE(ord("a")))
                plus = builder.call(strchr, (mode, INT(ord("+"))))
                plain = builder.icmp_unsigned("==", plus, ir.Constant(PTR, None))
                write_only = builder.and_(builder.or_(write, append), plain)
            else:
                access = builder.and_(mode, INT(os.O_ACCMODE))
                write_only = builder.icmp_unsigned("==", access, INT(os.O_WRONLY))
            builder.cbranch(write_only, allowed, inspect)
            builder.position_at_end(inspect)
        builder.cbranch(builder.call(check, (directory, path)), blocked, allowed)
        builder.position_at_end(blocked)
        builder.store(INT(13), builder.call(get_errno, ()))  # EACCES
        builder.ret(INT(-1) if result_type == INT else ir.Constant(PTR, None))
        builder.position_at_end(allowed)
        address = builder.call(lookup, (next_library, symbol))
        callee = builder.inttoptr(builder.ptrtoint(address, LONG), wrapper.function_type.as_pointer())
        builder.ret(builder.call(callee, wrapper.args))

    # The optional open mode is forwarded in the normal integer argument slot.
    for suffix in ("", "64"):
        wrap("open" + suffix, INT, (PTR, INT, INT))
        wrap("openat" + suffix, INT, (INT, PTR, INT, INT), at=True)
        wrap("__open" + suffix + "_2", INT, (PTR, INT))
        wrap("__openat" + suffix + "_2", INT, (INT, PTR, INT), at=True)
        wrap("fopen" + suffix, PTR, (PTR, PTR), stream=True)
        wrap("freopen" + suffix, PTR, (PTR, PTR, PTR), stream=True)
    wrap("opendir", PTR, (PTR,))
    return module


def build(directory):
    llvm.initialize_native_target()
    llvm.initialize_native_asmprinter()
    module = llvm.parse_assembly(str(create_filter()))
    module.verify()
    target = llvm.Target.from_default_triple().create_target_machine(reloc="pic")
    obj = directory / "file_access.o"
    obj.write_bytes(target.emit_object(module))
    linker = Path(sys.prefix) / "compiler_compat" / "ld"
    subprocess.run([linker, "-shared", "-z", "noexecstack", obj,
                    "-o", directory / "file_access.so"], check=True)


if __name__ == "__main__":
    build(Path(sys.argv[1]))

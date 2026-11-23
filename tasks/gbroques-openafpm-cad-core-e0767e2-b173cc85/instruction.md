Revise the OpenAFPM wind turbine so a builder can specify the clearance between the tail vane and the blade sweep at full furl. The linked baseline at /app/inputs/openafpm_cad_core has its production assembly at /app/inputs/openafpm_cad_core/documents/WindTurbine.FCStd. Deliver the complete revised Python and FreeCAD project at /app/out/openafpm_cad_core, with documents/WindTurbine.FCStd saved at rest for T Shape with both inputs at their defaults.

Supply one high-end stop for a family of interchangeable tail booms. Start from T Shape, H Shape or Star Shape in /app/inputs/openafpm_cad_core/get_default_parameters.py. Let B₀ be its preset BoomLength in mm. The builder may set the existing furling input BoomLength to any integer B in [0.60 B₀, 1.40 B₀], default B₀. Keep the other preset inputs. VaneClearanceFraction is an optional dimensionless furling input accepting any real q in [0.10, 0.15], default 0.10. For each (B, q), fit one maximum furl angle and one transverse stop width for all real service lengths b between B₀ and B, inclusive. Supply the assembly and fabrication tables for cut length B.

For each preset, let R be half its RotorDiameter in mm. The blade sweep is a filled, zero-thickness circular disk of radius R, centred at the mounted blade assembly back disk's native origin in its local z=0 plane. From the production rest pose, positive furl is rigid right-handed rotation of the moving tail about the mounted outer hinge's local +Z axis through its origin. Separations are minimum Euclidean distances between placed BReps in WindTurbine coordinates, measured with FreeCAD Part.

First blocking contact is the smallest furl angle in 0–130° at which the finished high-end stop is within 0.000001 mm of the yaw pipe. At that angle, resolved to 0.000001°, contact must be on the pipe's outer cylindrical wall. At 0.1° beyond contact, the stop and pipe must share positive volume. For clearance at service length b, translate the entire finished vane from its installed rest position at B by b−B mm along the mounted boom pipe's local +Z (toward the vane). It furls from there about the same hinge. At first contact the least vane-to-disk clearance over the whole service interval must be within 0.05 mm of qR. Every service length must retain at least qR − 0.05 mm throughout the preceding stroke. The native unrounded travel limit must agree with the contact angle within 0.01°. The motion export must report that native limit rounded to two decimal places in degrees.

The existing entry points must work from a relocated copy with its own linked documents. Assembly and table calls run separately, each in a fresh /opt/conda/bin/python process from the package's parent directory. Here p is the preset name. Omit either assignment to args[1] to use that input's default.

    from openafpm_cad_core.get_default_parameters import get_default_parameters
    from openafpm_cad_core.load import load_turbine
    from openafpm_cad_core.furl_transform import get_furl_transform
    from openafpm_cad_core.dimension_tables import load_dimension_tables
    cfg = get_default_parameters(p)
    args = [cfg['magnafpm'], dict(cfg['furling']), cfg['user']]
    args[1]['BoomLength'] = B
    args[1]['VaneClearanceFraction'] = q
    w, m = load_turbine(*args)
    motion = get_furl_transform(w, m)

Call load_dimension_tables(*args) with the same inputs. The saved T assembly must reopen and recompute directly in FreeCAD 1.1.0. With B = B₀, setting the cell aliased VaneClearanceFraction on Spreadsheet in documents/Master_of_Puppets.FCStd to any admissible q and recomputing the linked documents must yield the specified build with every object valid and none left touched.

For preservation, the comparator is the supplied project at the same preset and cut length B, recomputed with the released unrounded angle and width substituted for HighEndStop.MaximumFurlAngle and HighEndStop.HighEndStopWidth. Installed physical occurrences must correspond one for one and retain their moving or fixed relationships during furl. Require valid closed solids, material boundaries within 0.05 mm in both directions, and datum orientations within 0.01°. The published motion transforms must agree with the assembled datums to the same linear and angular tolerances. The native and tabulated fabrication dimensions must match the comparator at the tables' existing displayed precision. Every build across the stated domain must meet these requirements as returned and after recomputation.

Write only your final deliverables under `/app/out`.

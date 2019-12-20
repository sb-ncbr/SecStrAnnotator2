﻿using System;
using System.Collections.Generic;
using System.Linq;
using Cif.Components;

namespace protein.SSEs
{
	public enum SSEType
	{
		MIXED_HELIX_TYPE,
		HELIX_G_TYPE,
		HELIX_H_TYPE,
		HELIX_I_TYPE,
		SHEET_TYPE,
		ISOLATED_BETA_BRIDGE_TYPE,
		TURN_C7_TYPE,
		WIGGLE_C7_TYPE,
		BULGE_CLASSIC_SHORT_SIDE_TYPE,
		BULGE_CLASSIC_LONG_SIDE_TYPE,
		BULGE_WIDE_SHORT_SIDE_TYPE,
		BULGE_WIDE_LONG_SIDE_TYPE,
		BULGE_ANTIPARALLEL33_SHORT_SIDE_TYPE, // in 2qad chain B ~ resi 15 // "short" side is the one donating protons
		BULGE_ANTIPARALLEL33_LONG_SIDE_TYPE, // "long" side is the one accepting protons
		BULGE_ANTIPARALLEL22_SHORT_SIDE_TYPE, // in 1gei ~ resi 13 // "short" side is the one donating protons
		BULGE_ANTIPARALLEL22_LONG_SIDE_TYPE, // "long" side is the one accepting protons
		BULGE_ANTIPARALLEL15_SHORT_SIDE_TYPE, // in 1gjm ~ resi 94
		BULGE_ANTIPARALLEL15_LONG_SIDE_TYPE, 
		BULGE_ANTIPARALLEL23_SHORT_SIDE_TYPE, // in 3dbg ~ resi 452
		BULGE_ANTIPARALLEL23_LONG_SIDE_TYPE, 
		BULGE_PARALLEL14_SHORT_SIDE_TYPE,
		BULGE_PARALLEL14_LONG_SIDE_TYPE,
		BULGE_PARALLEL32_SHORT_SIDE_TYPE, // in 3ruk ~ resi 38
		BULGE_PARALLEL32_LONG_SIDE_TYPE,
		BULGE_PARALLEL13_SHORT_SIDE_TYPE, // in 3dax ~ resi 35
		BULGE_PARALLEL13_LONG_SIDE_TYPE,
		BULGE_PARALLEL33_SHORT_SIDE_TYPE, // in 3v8d ~ resi 69
		BULGE_PARALLEL33_LONG_SIDE_TYPE,
	}
}
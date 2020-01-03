using System;
using System.Collections.Generic;
using System.Linq;

using Cif.Components;
using protein.Libraries;

namespace protein.Sses
{
    static class LibSseTypes
    {
        public static readonly HashSet<SseType> ALL_HELIX_TYPES = new HashSet<SseType>{
            SseType.HELIX_H_TYPE,
            SseType.HELIX_G_TYPE,
            SseType.HELIX_I_TYPE,
            SseType.MIXED_HELIX_TYPE
        };

        public static readonly HashSet<SseType> ALL_SHEET_TYPES = new HashSet<SseType>{
            SseType.SHEET_TYPE,
            SseType.ISOLATED_BETA_BRIDGE_TYPE
        };

        private static readonly Dictionary<SseType, string> sseTypeToString = new Dictionary<SseType, string>{
            { SseType.NOT_FOUND_TYPE, "X" },

            { SseType.MIXED_HELIX_TYPE, "h" },
            { SseType.HELIX_G_TYPE, "G" },
            { SseType.HELIX_H_TYPE, "H" },
            { SseType.HELIX_I_TYPE, "I" },

            { SseType.SHEET_TYPE, "E" },
            { SseType.ISOLATED_BETA_BRIDGE_TYPE, "B" },

            { SseType.TURN_C7_TYPE, "C" },
            { SseType.WIGGLE_C7_TYPE, "W" },

            { SseType.BULGE_ANTIPARALLEL_UNSPECIFIED, "bA" },
            { SseType.BULGE_PARALLEL_UNSPECIFIED, "bP" },

            { SseType.BULGE_ANTIPARALLEL_CLASSIC_SHORT_SIDE, "n" },
            { SseType.BULGE_ANTIPARALLEL_CLASSIC_LONG_SIDE, "N" },
            { SseType.BULGE_ANTIPARALLEL_WIDE_SHORT_SIDE, "m" },
            { SseType.BULGE_ANTIPARALEL_WIDE_LONG_SIDE, "M" },
            { SseType.BULGE_ANTIPARALLEL33_SHORT_SIDE, "u" }, // in 2qad chain B ~ resi 15 // "short" side is the one donating protons
            { SseType.BULGE_ANTIPARALLEL33_LONG_SIDE, "U" }, // "long" side is the one accepting protons
            { SseType.BULGE_ANTIPARALLEL22_DONOR_SIDE, "t" }, // in 1gei ~ resi 13 // "short" side is the one donating protons
            { SseType.BULGE_ANTIPARALLEL22_ACCEPTOR_SIDE, "T" }, // "long" side is the one accepting protons
            { SseType.BULGE_ANTIPARALLEL15_SHORT_SIDE, "s" }, // in 1gjm ~ resi 94 -- very rare (8 in PDB)!
            { SseType.BULGE_ANTIPARALLEL15_LONG_SIDE, "S" },
            { SseType.BULGE_ANTIPARALLEL23_SHORT_SIDE, "o" }, // in 3dbg ~ resi 452
            { SseType.BULGE_ANTIPARALLEL23_LONG_SIDE, "O" },

            { SseType.BULGE_PARALLEL14_SHORT_SIDE, "p" },
            { SseType.BULGE_PARALLEL14_LONG_SIDE, "P" },
            { SseType.BULGE_PARALLEL32_SHORT_SIDE, "q" }, // in 3ruk ~ resi 38
            { SseType.BULGE_PARALLEL32_LONG_SIDE, "Q" },
            { SseType.BULGE_PARALLEL13_SHORT_SIDE, "r" }, // in 3dax ~ resi 35
            { SseType.BULGE_PARALLEL13_LONG_SIDE, "R" },
            { SseType.BULGE_PARALLEL33_SHORT_SIDE, "l" }, // in 3v8d ~ resi 69
            { SseType.BULGE_PARALLEL33_LONG_SIDE, "L" },
        };
        private static readonly Dictionary<string, SseType> stringToSseType = sseTypeToString.ToDictionary(kv => kv.Value, kv => kv.Key);



        public static string AsString(this SseType type, bool full = false) => full ? type.ToString() : sseTypeToString[type];

        public static SseType Type(string typeString) => stringToSseType[typeString];
        

        public static bool IsHelix(this SseType type) => ALL_HELIX_TYPES.Contains(type);

        public static bool IsSheet(this SseType type) => ALL_SHEET_TYPES.Contains(type);
    }
}

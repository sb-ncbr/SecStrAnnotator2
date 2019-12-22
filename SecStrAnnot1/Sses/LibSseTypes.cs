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

            { SseType.BULGE_CLASSIC_SHORT_SIDE_TYPE, "n" },
            { SseType.BULGE_CLASSIC_LONG_SIDE_TYPE, "N" },
            { SseType.BULGE_WIDE_SHORT_SIDE_TYPE, "m" },
            { SseType.BULGE_WIDE_LONG_SIDE_TYPE, "M" },
            { SseType.BULGE_ANTIPARALLEL33_SHORT_SIDE_TYPE, "u" }, // in 2qad chain B ~ resi 15 // "short" side is the one donating protons
            { SseType.BULGE_ANTIPARALLEL33_LONG_SIDE_TYPE, "U" }, // "long" side is the one accepting protons
            { SseType.BULGE_ANTIPARALLEL22_SHORT_SIDE_TYPE, "t" }, // in 1gei ~ resi 13 // "short" side is the one donating protons
            { SseType.BULGE_ANTIPARALLEL22_LONG_SIDE_TYPE, "T" }, // "long" side is the one accepting protons
            { SseType.BULGE_ANTIPARALLEL15_SHORT_SIDE_TYPE, "s" }, // in 1gjm ~ resi 94
            { SseType.BULGE_ANTIPARALLEL15_LONG_SIDE_TYPE, "S" },
            { SseType.BULGE_ANTIPARALLEL23_SHORT_SIDE_TYPE, "o" }, // in 3dbg ~ resi 452
            { SseType.BULGE_ANTIPARALLEL23_LONG_SIDE_TYPE, "O" },

            { SseType.BULGE_PARALLEL14_SHORT_SIDE_TYPE, "p" },
            { SseType.BULGE_PARALLEL14_LONG_SIDE_TYPE, "P" },
            { SseType.BULGE_PARALLEL32_SHORT_SIDE_TYPE, "q" }, // in 3ruk ~ resi 38
            { SseType.BULGE_PARALLEL32_LONG_SIDE_TYPE, "Q" },
            { SseType.BULGE_PARALLEL13_SHORT_SIDE_TYPE, "r" }, // in 3dax ~ resi 35
            { SseType.BULGE_PARALLEL13_LONG_SIDE_TYPE, "R" },
            { SseType.BULGE_PARALLEL33_SHORT_SIDE_TYPE, "l" }, // in 3v8d ~ resi 69
            { SseType.BULGE_PARALLEL33_LONG_SIDE_TYPE, "L" },
        };
        private static readonly Dictionary<string, SseType> stringToSseType = sseTypeToString.ToDictionary(kv => kv.Value, kv => kv.Key);



        public static string AsString(this SseType type, bool full = false) => full ? type.ToString() : sseTypeToString[type];

        public static SseType Type(string typeString) => stringToSseType[typeString];
        

        public static bool IsHelix(this SseType type) => ALL_HELIX_TYPES.Contains(type);

        public static bool IsSheet(this SseType type) => ALL_SHEET_TYPES.Contains(type);
    }
}

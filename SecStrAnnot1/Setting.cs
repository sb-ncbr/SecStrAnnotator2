#define DEVEL // DEVEL==true for development versions, DEVEL==false for release versions

using System;
using System.Collections.Generic;
using System.Linq;

using protein.Sses;
using protein.Libraries;

namespace protein
{
    static class Setting
    {
        public const string NAME = "SecStrAnnotator";
        public const string VERSION_NUMBER = "2.2";
#if DEVEL
        public static string VERSION = $"{VERSION_NUMBER}.{Lib.BuildVersion.Build}.{Lib.BuildVersion.Revision} [{Lib.BuildTime:u}]";
#else
		public static string VERSION = $"{VERSION_NUMBER} [{Lib.BuildTime:yyyy-MM-dd}]";
#endif

        public static string[] CommandLineArguments;

        public enum AlignMethod { None, Align, Super, Cealign };
        public static Dictionary<AlignMethod, string> alignMethodNames = new Dictionary<AlignMethod, string>
            { { AlignMethod.None,"none" }, { AlignMethod.Align,"align" }, { AlignMethod.Super,"super"}, { AlignMethod.Cealign,"cealign" }};
        public const AlignMethod DEFAULT_ALIGN_METHOD = AlignMethod.Cealign;

        public enum SecStrMethod { File, Dssp, Hbond1, Hbond2, Geom, GeomDssp, GeomHbond1, GeomHbond2 };

        public static Dictionary<SecStrMethod, string> secStrMethodNames = new Dictionary<SecStrMethod, string>{
            { SecStrMethod.File,"file" },
            { SecStrMethod.Dssp,"dssp" },
            { SecStrMethod.Hbond1,"hbond1" },
            { SecStrMethod.Hbond2,"hbond2" },
            { SecStrMethod.Geom,"geom" },
            { SecStrMethod.GeomDssp,"geom-dssp" },
            { SecStrMethod.GeomHbond1,"geom-hbond1" },
            { SecStrMethod.GeomHbond2,"geom-hbond2" } };
        public const SecStrMethod DEFAULT_SEC_STR_METHOD = SecStrMethod.GeomHbond2;

        public static SseType[] DEFAULT_ACCEPTED_SSE_TYPES = LibSseTypes.ALL_HELIX_TYPES.Concat(LibSseTypes.ALL_SHEET_TYPES).ToArray();

        public enum SelectionMethod { DynProg, BB, MOM, Combined, None };
        public static Dictionary<SelectionMethod, string> selectionMethodNames = new Dictionary<SelectionMethod, string>
        { { SelectionMethod.DynProg,"dp" }, { SelectionMethod.BB,"bb" }, { SelectionMethod.MOM,"mom" }, { SelectionMethod.Combined,"combined" }, { SelectionMethod.None,"none" } };
        public const SelectionMethod DEFAULT_SELECTION_METHOD = SelectionMethod.MOM;

        public enum MetricMethod { No3, No7, No8 };
        public static Dictionary<MetricMethod, string> metricMethodNames = new Dictionary<MetricMethod, string>
        { { MetricMethod.No3,"3" }, { MetricMethod.No7,"7" }, { MetricMethod.No8,"8" } };
        public const MetricMethod DEFAULT_METRIC_METHOD = MetricMethod.No8;


        // penalty for not matching template SSE S1 with query SSE S2 = pen0 + pen1*L1 + pen2*L2, where L1 (L2) is length of S1 (S2) in Angstroms.
        public static double[] DEFAULT_MAXMETRIC = new double[] { 30, 0, 0 };

        public const int DEFAULT_MAX_GAP_FOR_SOFT_MATCHING = 0;
        public const int DEFAULT_EXTRA_RESIDUES_IN_SEQUENCE = 0;

        public const double DEFAULT_RMSD_LIMIT = 1.00;

        public const double DEFAULT_H_BOND_ENERGY_LIMIT = -0.5;

        public const double STR_ALIGNMENT_SCALING_DISTANCE = 20.0; //Angstrom

        public const string DEFAULT_CHAIN_ID = "A";
        public const string DEFAULT_DOMAIN_RANGES_STRING = ":";

        public const int JSON_OUTPUT_MAX_INDENT_LEVEL = 3;

        public const bool LABEL_DETECTED_SSES_AS_NULL = false;

        public const string PDB_FILE_EXT = ".cif";
        public const string ALIGNED_PDB_FILE_EXT = "-aligned.cif";
        public const string ALIGNMENT_FILE_EXT = "-alignment.json";
        public const string RENUMBERED_PDB_FILE_EXT = "-renumbered.cif";
        public const string TEMPLATE_ANNOTATION_FILE_EXT = "-template.sses.json";
        public const string ANNOTATION_FILE_EXT = "-annotated.sses.json";
        public const string ANNOTATION_WITH_SEQUENCES_FILE_EXT = "-annotated_with_sequences.sses.json";
        public const string DSSP_OUTPUT_FILE_EXT = ".dssp";
        public const string LINE_SEGMENTS_FILE_EXT = "-line_segments.cif";
        public const string INPUT_SSES_FILE_EXT = ".sses.json";
        public const string DETECTED_SSES_FILE_EXT = "-detected.sses.json";
        public const string JOINED_SSES_FILE_EXT = "-joined.sses.json";
        public const string RMSDS_FILE_EXT = "-rmsds.tsv";
        public const string LABEL2AUTH_FILE_EXT = "-label2auth.tsv";
        public const string SEQUENCE_FILE_EXT = "-sequence.json";

        public const string CONFIG_FILE = "SecStrAnnotator_config.json";

        public const bool FILTER_OUTPUT_BY_LABEL = false;
        public static string[] OUTPUT_ONLY_THESE_LABELS = new string[] { "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "F'", "G'", "J'", "K'", "1a", "1b", "1c", "1d", "2a", "2b" };

        public static bool IgnoreInsertions { get; set; }
        public static bool IgnoreInsertionsWarningThrown { get; set; }

        // joiningTypeCombining(X,Y) should determine the type of SSE resulting from joining SSE of type X with SSE of type Y. 
        // Return value null indicates that SSEs of these 2 types cannot be joined. 
        public static Func<SseType, SseType, SseType?> JoiningTypeCombining =>
            (x, y) =>
				(x == y) ? x
				: (x.IsHelix() && y.IsHelix()) ? SseType.MIXED_HELIX_TYPE
				: (x.IsSheet() && y.IsSheet()) ? SseType.SHEET_TYPE
				: (SseType?)null;

        public static string Directory { get; set; }

    }
}
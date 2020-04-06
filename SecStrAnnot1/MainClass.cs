using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Resources;
using System.Threading;
using System.Globalization;

using SecStrAnnotator2.Utils;
using Cif.Components;
using protein.Libraries;
using protein.Sses;
using protein.SecStrAssigning;
using protein.Annotating;
using protein.Json;

namespace protein
{
    class MainClass
    {
        public static int Main_SecStrAnnot1(string[] args)
        {
            #region Declarations.
            Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;
            Thread.CurrentThread.CurrentUICulture = CultureInfo.InvariantCulture;

            List<(String, TimeSpan)> times = new List<(String, TimeSpan)>();
            DateTime t0 = DateTime.Now;
            DateTime stamp = DateTime.Now;

            String templateID = null;
            String queryID = null;

            bool onlyDetect = false;

            string templateChainID_ = Setting.DEFAULT_CHAIN_ID;
            string queryChainID_ = Setting.DEFAULT_CHAIN_ID;

            List<(int, int)> templateDomainRanges = null;
            List<(int, int)> queryDomainRanges = null;

            Setting.AlignMethod alignMethod = Setting.DEFAULT_ALIGN_METHOD;
            bool tryToReuseAlignment = false;

            Setting.SecStrMethod secStrMethod = Setting.DEFAULT_SEC_STR_METHOD;
            bool joinHelices = false;
            bool forceCalculateVectors = false;

            Setting.SelectionMethod selectionMethod = Setting.DEFAULT_SELECTION_METHOD;
            bool alternativeMatching = false;
            bool softMatching = false;
            int maxGapForSoftMatching = Setting.DEFAULT_MAX_GAP_FOR_SOFT_MATCHING;
            int extraResiduesInSequence = Setting.DEFAULT_EXTRA_RESIDUES_IN_SEQUENCE;

            string correctionsFile = null;

            bool createPymolSession = false;

            SseType[] acceptedSSETypes = Setting.DEFAULT_ACCEPTED_SSE_TYPES;
            //LibAnnotation.JoiningParameters joiningParameters = new LibAnnotation.JoiningParameters (LibAnnotation.ALL_HELIX_TYPES, 10, 1, 10, 2, 30);
            LibAnnotation.JoiningParameters joiningParameters = new LibAnnotation.JoiningParameters(LibSseTypes.ALL_HELIX_TYPES.ToArray(), 10, 1, 0, 5, 60); //try stagger penalty 10 or 5

            double rmsdLimit = Setting.DEFAULT_RMSD_LIMIT;

            // annotationTypeFitting(T,Q) should determine whether query SSE of type Q can be mapped to template SSE of type T.
            Func<SseType, SseType, bool> annotationTypeFitting =
                (t, q) => t.IsHelix() && q.IsHelix() || t.IsSheet() && q.IsSheet();

            Setting.MetricMethod metricMethod = Setting.DEFAULT_METRIC_METHOD;
            double[] maxmetric = Setting.DEFAULT_MAXMETRIC;
            Func<SseInSpace, double> skipTemplatePenalty = (sse => maxmetric[0] + maxmetric[1] * (sse.EndPoint - sse.StartPoint).Size);
            Func<SseInSpace, double> skipCandidatePenalty = (sse => maxmetric[2] * (sse.EndPoint - sse.StartPoint).Size);

            double? momTimeoutSeconds = null;

            bool printLabel2AuthTable = false;

            #endregion


            #region Processing options and arguments.
            Options options = new Options();
            options.GlobalHelp = Setting.NAME + " " + Setting.VERSION;

            options.AddArgument(new Argument("DIRECTORY")
                .AddHelp("Directory for all input and output files.")
            );
            options.AddArgument(new Argument("TEMPLATE")
                .AddHelp("Template domain specification in format PDB[,CHAIN,[RANGES]]")
            );
            options.AddArgument(new Argument("QUERY")
                .AddHelp("Query domain specification in format PDB[,CHAIN,[RANGES]]")
                .AddHelp("")
                .AddHelp("Domain specification examples: 1tqn or 1og2,B or 1tqn,A,28: or 1h9r,A,123:182,255:261")
                .AddHelp("Default CHAIN is A, default RANGES is : (whole chain).")
            );

            double _;

            options.AddOption(Option.DictionaryChoiceOption(new string[] { "-a", "--align" }, v => { alignMethod = v; }, Setting.alignMethodNames)
                .AddParameter("METHOD")
                .AddHelp("Specify structure alignment method.")
                .AddHelp("METHOD is one of " + Setting.alignMethodNames.Values.EnumerateWithCommas() + "; default: " + Setting.alignMethodNames[Setting.DEFAULT_ALIGN_METHOD])
                .AddHelp("    none:    use the query structure as it is, without aligning")
#if DEVEL
				.AddHelp("    simple:  a naive alignment using only 3 points in heme cofactor")
#endif
                .AddHelp("    align:   run PyMOL's command align")
                .AddHelp("    super:   run PyMOL's command super")
                .AddHelp("    cealign: run PyMOL's command cealign")
            );
            options.AddOption(Option.DictionaryChoiceOption(new string[] { "-d", "--ssa" }, v => { secStrMethod = v; }, Setting.secStrMethodNames)
                .AddParameter("METHOD")
                .AddHelp("Specify method of secondary structure assignment (SSA).")
                .AddHelp("METHOD is one of " + Setting.secStrMethodNames.Values.EnumerateWithCommas() + "; default: " + Setting.secStrMethodNames[Setting.DEFAULT_SEC_STR_METHOD])
                .AddHelp("    file:       read from file <QUERY_ID>" + Setting.INPUT_SSES_FILE_EXT)
                .AddHelp("    dssp:       run DSSP")
                .AddHelp("    hbond1:      use built-in DSSP-like algorithm (old version)")
                .AddHelp("    hbond2:      use built-in DSSP-like algorithm")
                .AddHelp("    geom:       use built-in geometry-based method")
                .AddHelp("    geom-dssp:  use geom for helices, dssp for sheets")
                .AddHelp("    geom-hbond1: use geom for helices, hbond1 for sheets")
                .AddHelp("    geom-hbond2: use geom for helices, hbond2 for sheets")
            );
            options.AddOption(Option.SwitchOption(new string[] { "-o", "--onlyssa" }, v => { onlyDetect = v; })
                .AddHelp("Run only secondary structure assignment (SSA), do not run annotation.")
                .AddHelp("Changes usage: dotnet " + System.AppDomain.CurrentDomain.FriendlyName + ".dll --onlyssa [OPTIONS] DIRECTORY QUERY")
            );
            options.AddOption(Option.DoubleOption(new string[] { "-l", "--limit" }, v => { rmsdLimit = v; })
                .AddParameter("LIMIT")
                .AddHelp("Specify RMSD limit for geometry-based secondary structure assignment, default: " + Setting.DEFAULT_RMSD_LIMIT.ToString("0.0#####"))
            );
            options.AddOption(Option.StringOption(new string[] { "-t", "--types" }, v => { acceptedSSETypes = v.Split(',').Select(str => LibSseTypes.Type(str)).ToArray(); })
                .AddConstraint(optArgs => optArgs[0].Split(',').All(type => type.Length == 1), "must be a comma-separated list of one-character SSE types")
                .AddConstraint(optArgs => optArgs[0].Split(',').All(type => LibSseTypes.IsType(type)), "contains other than allowed types")
                .AddParameter("TYPES")
                .AddHelp("Specify the allowed types of secondary structures.")
                .AddHelp("TYPES is a list of comma-separated (without space) types; default: " + Setting.DEFAULT_ACCEPTED_SSE_TYPES.Select(type => type.AsString()).EnumerateWithSeparators(","))
                .AddHelp("The types are denoted by DSSP convention (H = alpha helix, G = 3_10 helix, I = pi helix, h = helix, E = strand with >2 H-bonds, B = strand with 2 H-bonds)")
            );
            options.AddOption(Option.DictionaryChoiceOption(new string[] { "-m", "--matching" }, v => { selectionMethod = v; }, Setting.selectionMethodNames)
                .AddParameter("METHOD")
                .AddHelp("Specify method of matching template to query SSEs.")
                .AddHelp("METHOD is one of " + Setting.selectionMethodNames.Values.EnumerateWithCommas() + "; default: " + Setting.selectionMethodNames[Setting.DEFAULT_SELECTION_METHOD])
                .AddHelp("    dp:      Dynamic Programming (fast, ignores beta-strand connectivity)")
                .AddHelp("    mom:     Mixed Ordered Matching, branch&bound max-weight-clique algorithm (follows beta-strand connectivity)")
#if DEVEL
				.AddHelp("    bb:      older branch&bound algorithm, max-weight-clique algorithm (slow, follows beta-strand connectivity)") 
				.AddHelp("    combined:runs dynamic programming first; if connectivity issues arise, runs branch&bound for sheets only, than reruns DP with constraints from BB (faster that pure BB, follows beta-strand connectivity, might give suboptimal solution)")
				.AddHelp("    none:    do not run selection algorithm, only compute metric matrix")
#endif
            );
            options.AddOption(Option.SwitchOption(new string[] { "-n", "--soft" }, v => { softMatching = v; })
#if DEVEL
				.AddHelp("Use soft matching in MOM algorithm (and BB), i.e. allow matching 2 query strands to 1 template strand, et vice versa.")
#else
                .AddHelp("Use soft matching in MOM algorithm, i.e. allow matching 2 query strands to 1 template strand, et vice versa.")
#endif
            );
            options.AddOption(Option.SwitchOption(new string[] { "-s", "--session" }, v => { createPymolSession = v; })
                .AddHelp("Create PyMOL session with results (.pse file).")
            );
            options.AddOption(Option.DictionaryChoiceOption(new string[] { "-M", "--metrictype" }, v => { metricMethod = v; }, Setting.metricMethodNames)
                .AddParameter("TYPE")
                .AddHelp("Specify metric for measuring difference between two SSEs.")
                .AddHelp("TYPE is one of " + Setting.metricMethodNames.Values.EnumerateWithCommas() + "; default: " + Setting.metricMethodNames[Setting.DEFAULT_METRIC_METHOD])
                .AddHelp("    3:  based on 3D coordinates (distance of start vectors + distance of end vectors)")
                .AddHelp("    7:  based on residue positions in structural alignment")
                .AddHelp("    8:  0.5 * (metric3 + metric7) + length_difference_penalty")
            );
            options.AddOption(Option.StringOption(new string[] { "-k", "--maxmetric" }, v => { maxmetric = v.Split(',').Select(Double.Parse).ToArray().Resized(3); })
                .AddConstraint(optArgs => optArgs[0].Split(',').All(k => Double.TryParse(k, out _)), "must be a float or a comma-separated list 3 floats")
                .AddParameter("K0[,K1,K2]")
                .AddHelp("Specify maximum allowed metric value K for matching template SSE S1 to query SSE S2.")
                .AddHelp("K = K0 + K1*L1 + K2*L2, where L1, L2 are lengths of S1, S2 in Angstroms; default: " + maxmetric.EnumerateWithSeparators(","))
            );
            options.AddOption(Option.DoubleOption(new string[] { "-f", "--fallback" }, v => { momTimeoutSeconds = v; })
                .AddParameter("LIMIT")
                .AddHelp("Specify time limit (in seconds) for MOM algorithm; in case of timeout falls back to DP algorithm")
            );
            options.AddOption(Option.SwitchOption(new string[] { "-c", "--calculatevectors" }, v => { forceCalculateVectors = v; })
                .AddHelp("Force calculation of start and end vector of SSEs, even if they are in the input files.")
            );
            options.AddOption(Option.SwitchOption(new string[] { "-v", "--verbose" }, v => { Lib.DoWriteDebug = v; })
                .AddHelp("Print debugging notes and write out additional files.")
            );
            options.AddOption(Option.SwitchOption(new string[] { "-L", "--label2auth" }, v => { printLabel2AuthTable = v; })
                .AddHelp("Create an additional file with a table for converting label_ to auth_ numbering.")
            );

#if DEVEL
            options.AddOption(Option.SwitchOption(new string[] { "-i", "--ignoreinsertions" }, v => { Setting.IgnoreInsertions = v; })
                .AddHelp("Ignore residues with insertion code, if such are present in the input file.")
                .AddHelp("WARNING: Loaded structure will not correspond fully to the input file!")
            );
			options.AddOption (Option.StringOption(new string[]{"-C", "--corrections"}, v => { correctionsFile = v; })
				.AddParameter ("CORR_FILE")
				.AddHelp("Override automatic annotations by annotations from file CORR_FILE.")
				.AddHelp("Required file format: tab-separated with columns [pdbId, label, chainId, start, end] (start=end=0 for NotFound).")
			);
			options.AddOption (Option.IntOption(new string[]{"-e", "--extraresidues"}, v => { extraResiduesInSequence = v; })
				.AddParameter ("N")
				.AddHelp("Specify the number of residues before and after each SSE that should be included in its sequence in output.")
			);
			options.AddOption (Option.IntOption(new string[]{"-g", "--maxgap"}, v => { maxGapForSoftMatching = v; })
				.AddParameter ("GAP")
				.AddHelp("Specify maximum gap for strand joining in soft matching, default: " + DEFAULT_MAX_GAP_FOR_SOFT_MATCHING.ToString())
			);
			options.AddOption (Option.SwitchOption(new string[]{"-j", "--join"}, v => { joinHelices = v; })
				.AddHelp("Refine detected SSEs by joining if they fulfil certain criteria.")
			);
			options.AddOption (Option.SwitchOption(new string[]{"-N", "--alternativematching"}, v => { alternativeMatching = v; })
				.AddHelp("Use alternative matching in BB algorithm, i.e. read alternatively mergeable template SSEs from template annotation.")
			);
			options.AddOption (Option.SwitchOption(new string[]{"-r", "--reusealignment"}, v => { tryToReuseAlignment = v; })
				.AddHelp("Reuse aligned query structure from <QUERY_ID>" + ALIGNED_PDB_FILE_EXT + ", if it exists.")
			);
#endif

            List<String> otherArgs;
            bool optionsOK = options.TryParse(args, out otherArgs);
            if (!optionsOK)
            {
                Environment.Exit(1);
            }

            string templateDomainString;
            string queryDomainString;
            if (onlyDetect)
            {
                // Execution with --onlyssa: SecStrAnnotator.exe DIR QUERY
                if (otherArgs.Count != 2)
                {
                    Options.PrintError("Exactly 2 arguments required when run with --onlyssa ({0} given: {1})", otherArgs.Count, otherArgs.EnumerateWithSeparators(" "));
                    Environment.Exit(1);
                }
                Setting.Directory = otherArgs[0];
                queryDomainString = otherArgs[1];
                try
                {
                    (queryID, queryChainID_, queryDomainRanges) = ParseDomainSpecification(queryDomainString);
                }
                catch (FormatException e)
                {
                    Options.PrintError("Invalid value of the 2nd argument \"" + queryDomainString + "\"\n"
                        + e.Message);
                    Environment.Exit(1);
                }
            }
            else
            {
                // Normal execution: SecStrAnnotator.exe DIR TEMPLATE TEMPLATE QUERY
                if (otherArgs.Count != 3)
                {
                    Options.PrintError("Exactly 3 arguments required ({0} given: {1})", otherArgs.Count, otherArgs.EnumerateWithSeparators(" "));
                    Environment.Exit(1);
                }
                Setting.Directory = otherArgs[0];
                templateDomainString = otherArgs[1];
                queryDomainString = otherArgs[2];
                try
                {
                    (templateID, templateChainID_, templateDomainRanges) = ParseDomainSpecification(templateDomainString, defaultChain: Setting.DEFAULT_CHAIN_ID);
                }
                catch (FormatException e)
                {
                    Options.PrintError("Invalid value of the 2nd argument \"" + templateDomainString + "\"\n"
                        + e.Message);
                    Environment.Exit(1);
                }
                try
                {
                    (queryID, queryChainID_, queryDomainRanges) = ParseDomainSpecification(queryDomainString, defaultChain: Setting.DEFAULT_CHAIN_ID);
                }
                catch (FormatException e)
                {
                    Options.PrintError("Invalid value of the 3rd argument \"" + queryDomainString + "\"\n"
                        + e.Message);
                    Environment.Exit(1);
                }
            }
            #endregion


            #region Reading configuration file.
            String configFile = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, Setting.CONFIG_FILE);
            Config config = new Config(configFile);
            #endregion


            #region Creating names for input and output files.
            String fileTemplatePDB = Path.Combine(Setting.Directory, templateID + Setting.PDB_FILE_EXT);
            String fileTemplateAnnotatedHelices = Path.Combine(Setting.Directory, templateID + Setting.TEMPLATE_ANNOTATION_FILE_EXT) + (Setting.JSON_INPUT ? ".json" : "");

            String fileQueryPDB = Path.Combine(Setting.Directory, queryID + Setting.PDB_FILE_EXT);
            String fileQueryAlignedPDB = Path.Combine(Setting.Directory, queryID + Setting.ALIGNED_PDB_FILE_EXT);
            String fileQueryAlignment = Path.Combine(Setting.Directory, queryID + Setting.ALIGNMENT_FILE_EXT);
            String fileQueryRenumberedPDB = Path.Combine(Setting.Directory, queryID + Setting.RENUMBERED_PDB_FILE_EXT);
            String fileQueryDSSP = Path.Combine(Setting.Directory, queryID + Setting.DSSP_OUTPUT_FILE_EXT);
            String fileQueryInputHelices = Path.Combine(Setting.Directory, queryID + Setting.INPUT_SSES_FILE_EXT) + (Setting.JSON_INPUT ? ".json" : "");
            String fileQueryDetectedHelices = Path.Combine(Setting.Directory, queryID + Setting.DETECTED_SSES_FILE_EXT) + (Setting.JSON_OUTPUT ? ".json" : "");
            String fileQueryJoinedHelices = Path.Combine(Setting.Directory, queryID + Setting.JOINED_SSES_FILE_EXT) + (Setting.JSON_OUTPUT ? ".json" : "");
            String fileQueryAnnotatedHelices = Path.Combine(Setting.Directory, queryID + Setting.ANNOTATION_FILE_EXT) + (Setting.JSON_OUTPUT ? ".json" : "");
            String fileQueryRmsds = Path.Combine(Setting.Directory, queryID + Setting.RMSDS_FILE_EXT);
            String fileQueryLabel2Auth = Path.Combine(Setting.Directory, queryID + Setting.LABEL2AUTH_FILE_EXT);

            times.Add(("Process options, arguments, create file names.", DateTime.Now.Subtract(stamp)));
            stamp = DateTime.Now;
            #endregion


            #region Reading proteins from PDB files.

            Protein tProtein; // template
            Protein qProtein; // query

            if (!onlyDetect)
            {
                // Lib.WriteInColor(ConsoleColor.Yellow, "Template protein:\n");
                tProtein = ReadProteinFromFile(fileTemplatePDB, templateChainID_, templateDomainRanges);
            }
            else
            {
                tProtein = new Protein(null as Cif.Tables.Model); // dummy protein, never will be used
            }

            // Lib.WriteInColor(ConsoleColor.Yellow, "Query protein:\n");
            if (tryToReuseAlignment && File.Exists(fileQueryAlignedPDB))
            {
                qProtein = ReadProteinFromFile(fileQueryAlignedPDB, queryChainID_, queryDomainRanges);
            }
            else
            {
                qProtein = ReadProteinFromFile(fileQueryPDB, queryChainID_, queryDomainRanges);
            }

            if (printLabel2AuthTable)
            {
                qProtein.SaveLabel2AuthTable(fileQueryLabel2Auth);
            }


            // string[] queryChainIDs = new string[] { queryChainID_ };
            // Dictionary<string, string[]> chainMapping = new Dictionary<string, string[]> { { templateChainID_, queryChainIDs } }; // mapping which chain in template corresponds to which chain(s) in the query protein
            // string[] allTemplateChainIDs = chainMapping.Keys.ToArray();
            // string[] allQueryChainIDs = chainMapping.Values.SelectMany(x => x).ToArray();
            // if (!onlyDetect)
            // {
            //     foreach (string c in chainMapping.Keys)
            //     {
            //         if (!tProtein.HasChain(c))
            //         {
            //             Lib.WriteError("Template protein does not contain chain {0}.", c);
            //             return -1;
            //         }
            //     }
            //     foreach (string c in allQueryChainIDs)
            //     {
            //         if (!qProtein.HasChain(c))
            //         {
            //             Lib.WriteError("Query protein does not contain chain {0}.", c);
            //             return -1;
            //         }
            //     }
            // }

            string[] queryChainIDs;// = new string[] { queryChainID_ };
                                   // Dictionary<string, string[]> chainMapping = new Dictionary<string, string[]> { { templateChainID_, queryChainIDs } }; // mapping which chain in template corresponds to which chain(s) in the query protein

            if (!onlyDetect)
            {
                if (!tProtein.HasChain(templateChainID_))
                {
                    Lib.WriteError("Template protein does not contain chain {0}.", templateChainID_);
                    return -1;
                }
            }

            if (queryChainID_ != null)
            {
                queryChainIDs = new string[] { queryChainID_ };
                // chainMapping = new Dictionary<string, string[]> { { templateChainID_, queryChainIDs } }; // mapping which chain in template corresponds to which chain(s) in the query protein
                foreach (string c in queryChainIDs)
                {
                    if (!qProtein.HasChain(c))
                    {
                        Lib.WriteError("Query protein does not contain chain {0}.", c);
                        return -1;
                    }
                }
            }
            else
            {
                queryChainIDs = qProtein.GetChains().Select(c => c.Id).ToArray();
                // chainMapping = new Dictionary<string, string[]> { { templateChainID_, queryChainIDs } }; // mapping which chain in template corresponds to which chain(s) in the query protein
            }
            // string[] allTemplateChainIDs = chainMapping.Keys.ToArray();
            // string[] allQueryChainIDs = chainMapping.Values.SelectMany(x => x).ToArray();



            #endregion


            #region Writing short info about the proteins to the console.

            times.Add(("Read files", DateTime.Now.Subtract(stamp)));
            stamp = DateTime.Now;
            #endregion


            #region Reading template annotation from file.
            SecStrAssignment templateSSA;
            if (!onlyDetect)
            {
                ISecStrAssigner templateSecStrAssigner =
                    Setting.JSON_INPUT ?
                    new FileSecStrAssigner_Json(fileTemplateAnnotatedHelices, templateID, new[] { templateChainID_ } /*allTemplateChainIDs*/)
                    : new FileSecStrAssigner(fileTemplateAnnotatedHelices, new[] { templateChainID_ } /*allTemplateChainIDs*/) as ISecStrAssigner;
                try
                {
                    templateSSA = templateSecStrAssigner.GetSecStrAssignment();
                }
                catch (SecStrAssignmentException)
                {
                    return -1;
                }
                catch (Exception)
                {
                    Lib.WriteErrorAndExit("Reading template annotation failed.");
                    return -1;
                }
                if (templateSSA.SSEs.Any(s => s.IsSheet && s.SheetId == null))
                {
                    foreach (Sse sse in templateSSA.SSEs.Where(s => s.IsSheet))
                    {
                        String num = String.Concat(sse.Label.TakeWhile(c => '0' <= c && c <= '9'));
                        if (num.Length == 0)
                        {
                            Lib.WriteWarning("Cannot determine sheet ID from label \"{0}\". Assigning sheet ID null.", sse.Label);
                            sse.SheetId = null;
                        }
                        else
                        {
                            sse.SheetId = Int32.Parse(num);
                        }
                    }
                }
            }
            else
            {
                templateSSA = null;
            }
            #endregion


            #region Secondary Structure Assignment of the query protein.

            ISecStrAssigner secStrAssigner;

            switch (secStrMethod)
            {
                case Setting.SecStrMethod.File:
                    secStrAssigner = Setting.JSON_INPUT ?
                        new FileSecStrAssigner_Json(fileQueryInputHelices, queryID, queryChainIDs)
                        : new FileSecStrAssigner(fileQueryInputHelices, queryChainIDs) as ISecStrAssigner;
                    break;
                case Setting.SecStrMethod.Dssp:
                    secStrAssigner = new DsspSecStrAssigner(qProtein, config.DsspExecutable, fileQueryRenumberedPDB, fileQueryDSSP, queryChainIDs, acceptedSSETypes);
                    break;
                case Setting.SecStrMethod.Geom:
                    secStrAssigner = new GeomSecStrAssigner(queryChainIDs.Select(c => qProtein.GetChain(c)), rmsdLimit);
                    break;
                case Setting.SecStrMethod.GeomDssp:
                    secStrAssigner = new GeomDsspSecStrAssigner(qProtein, queryChainIDs, rmsdLimit, config.DsspExecutable, fileQueryRenumberedPDB, fileQueryDSSP, acceptedSSETypes);
                    break;
                case Setting.SecStrMethod.Hbond1:
                    secStrAssigner = new HBondSecStrAssigner(qProtein, Setting.DEFAULT_H_BOND_ENERGY_LIMIT);
                    break;
                case Setting.SecStrMethod.GeomHbond1:
                    secStrAssigner = new GeomHbondSecStrAssigner(qProtein, rmsdLimit, Setting.DEFAULT_H_BOND_ENERGY_LIMIT);
                    break;
                case Setting.SecStrMethod.Hbond2:
                    secStrAssigner = new HBondSecStrAssigner2(qProtein, Setting.DEFAULT_H_BOND_ENERGY_LIMIT);
                    break;
                case Setting.SecStrMethod.GeomHbond2:
                    secStrAssigner = new GeomHbondSecStrAssigner2(qProtein, rmsdLimit, Setting.DEFAULT_H_BOND_ENERGY_LIMIT);
                    break;
                default:
                    throw new Exception("Unknown secondary structure detection method.");
            }

            if (joinHelices)
            {
                secStrAssigner = new JoiningSecStrAssigner(secStrAssigner, qProtein, rmsdLimit, Setting.JoiningTypeCombining);
            }

            secStrAssigner = new FilteringSecStrAssigner(secStrAssigner, acceptedSSETypes, queryChainIDs);
            if (secStrMethod != Setting.SecStrMethod.File)
            {
                secStrAssigner = new RelabellingSecStrAssigner(secStrAssigner, null, Setting.LABEL_DETECTED_SSES_AS_NULL);
            }
            secStrAssigner = new AuthFieldsAddingSecStrAssigner(secStrAssigner, qProtein);
            if (Setting.JSON_OUTPUT)
                secStrAssigner = new OutputtingSecStrAssigner_Json(secStrAssigner, fileQueryDetectedHelices, queryID);
            else
                secStrAssigner = new OutputtingSecStrAssigner(secStrAssigner, fileQueryDetectedHelices);

            SecStrAssignment querySSA;
            try
            {
                querySSA = secStrAssigner.GetSecStrAssignment();
            }
            catch (SecStrAssignmentException e)
            {
                Lib.WriteError(e.Message);
                return -1;
            }

            times.Add(("Secondary structure assignment", DateTime.Now.Subtract(stamp)));
            stamp = DateTime.Now;

            #endregion


            if (onlyDetect)
            {
                if (createPymolSession)
                {
                    if (!Setting.JSON_OUTPUT)
                    {
                        throw new NotImplementedException("Creating PyMOL session is implemented only for JSON output.");
                    }
                    Lib.WriteInColor(ConsoleColor.Yellow, "Running PyMOL:\n");
                    if (!Lib.RunPyMOLScriptWithCommandLineArguments(config.PymolExecutable, config.PymolScriptSession, new string[]{
                        "cif",
                        Setting.Directory,
                        queryDomainString,
                        "--detected",
                        "--hbonds",
                        }))
                        return -1;
                }
                Lib.WriteInColor(ConsoleColor.Yellow, "Detected SSEs: {0}\n", querySSA.SSEs.Count);
                if (Lib.DoWriteDebug)
                {
                    var qSSEsInSpace = LibAnnotation.SSEsAsLineSegments_GeomVersion(qProtein, querySSA.SSEs);
                    LibAnnotation.WriteAnnotationFile_Json(
                        fileQueryDetectedHelices,
                        queryID,
                        qSSEsInSpace,
                        null, //extras
                        querySSA.Connectivity,
                        querySSA.HBonds,
                        null,
                        null);
                    PrintTimes(times, t0);
                }
                return 0;
            }

            List<SseInSpace> annotQHelicesInSpace_AllChains = new List<SseInSpace>();
            List<(int, int, int)> annotQConnectivity_AllChains = new List<(int, int, int)>();
            List<double> metricList_AllChains = new List<double>();
            List<double> suspiciousnessList_AllChains = new List<double>();
            List<List<double>> rmsdLists_AllChains = new List<List<double>>();
            List<string> annotQSseSequences_AllChains = new List<string>();
            List<double> metric3List_AllChains = new List<double>();
            List<double> metric7List_AllChains = new List<double>();

            List<SseInSpace> detQSsesInSpace_AllChains = new List<SseInSpace>();
            List<String> detQSseSequences_AllChains = new List<string>();

            JsonValue rotationMatrix = null;

            #region Superimposition and annotation for each chain.
            // foreach (string templateChainID in allTemplateChainIDs)
            // {
            foreach (string queryChainID in queryChainIDs)
            {
                #region Superimposing p over t.
                if (tryToReuseAlignment && File.Exists(fileQueryAlignedPDB))
                {
                    // The query protein has already been read from aligned PDB file.
                }
                else
                {
                    if (alignMethod == Setting.AlignMethod.None)
                    {
                        // qProtein.Save (fileQueryAlignedPDB);
                        File.Copy(fileQueryPDB, fileQueryAlignedPDB, true);
                    }
                    else if (Setting.alignMethodNames.ContainsKey(alignMethod))
                    {
                        #region Superimpose by a given PyMOL's command.
                        Lib.WriteInColor(ConsoleColor.Yellow, "Running PyMOL to align proteins:\n");
                        if (!Lib.RunPyMOLScriptWithCommandLineArguments(config.PymolExecutable, config.PymolScriptAlign, new string[] {
                                "cif",
                                Setting.alignMethodNames [alignMethod],
                                Setting.Directory,
                                templateID,
                                templateChainID_.ToString (),
                                FormatRanges(templateDomainRanges),
                                queryID,
                                queryChainID.ToString (),
                                FormatRanges(queryDomainRanges)
                            }))
                            return -1;
                        qProtein = ReadProteinFromFile(fileQueryAlignedPDB, queryChainID_, queryDomainRanges).KeepOnlyNormalResidues(true);
                        rotationMatrix = LibAnnotation.ReadRotationMatrixFromAlignmentFile(fileQueryAlignment);
                        // Console.WriteLine(rotationMatrix);
                        #endregion
                    }
                    else
                    {
                        throw new Exception(alignMethod + " has no assigned name in alignMethodNames.");
                    }
                }
                times.Add(("Alignment", DateTime.Now.Subtract(stamp)));
                stamp = DateTime.Now;
                #endregion

                Lib.Shuffler shuffler;
                List<Sse> tSSEs = templateSSA.SSEs.WhereAndGetShuffler(sse => sse.ChainID == templateChainID_, out shuffler).ToList();
                List<(int, int, int)> tConnectivity = shuffler.UpdateIndices(templateSSA.Connectivity).ToList();
                List<Sse> qSSEs = querySSA.SSEs.WhereAndGetShuffler(sse => sse.ChainID == queryChainID, out shuffler).ToList();
                List<(int, int, int)> qConnectivity = shuffler.UpdateIndices(querySSA.Connectivity).ToList();
                List<SseInSpace> tSSEsInSpace;
                List<SseInSpace> qSSEsInSpace;


                #region Calculate helix-fit and sheet-fit RMSDs for whole chain and output them to a file.
                if (Lib.DoWriteDebug)
                {
                    const int UNIT_LENGTH = 4;
                    TextWriter wRMSD = new StreamWriter(fileQueryRmsds);
                    wRMSD.WriteLine(LibAnnotation.COMMENT_SIGN_WRITE + "RMSDs from fitting helix and sheet shape to the alpha-carbons. Value for resi=x is calculated from residues x, x+1, x+2, x+3.");
                    wRMSD.WriteLine(LibAnnotation.COMMENT_SIGN_WRITE + "resi\tRMSD_vs_H\tRMSD_vs_E");
                    List<Residue> residues = qProtein.GetChain(queryChainID).GetResidues().ToList();
                    for (int i = 0; i <= residues.Count - UNIT_LENGTH; i++)
                    {
                        double rmsdH;
                        double rmsdE;
                        // double rmsdH5;
                        // double rmsdH6;
                        LibAnnotation.CheckGeometryOf1Unit(residues.GetRange(i, UNIT_LENGTH), SseType.HELIX_H_TYPE, rmsdLimit, out rmsdH);
                        LibAnnotation.CheckGeometryOf1Unit(residues.GetRange(i, UNIT_LENGTH), SseType.SHEET_TYPE, rmsdLimit, out rmsdE);
                        // LibAnnotation.CheckGeometryOf1Unit(residues.GetRange(i, 5), '5', rmsdLimit, out rmsdH5);
                        // LibAnnotation.CheckGeometryOf1Unit(residues.GetRange(i, 6), '6', rmsdLimit, out rmsdH6);
                        wRMSD.WriteLine("{0}\t{1}\t{2}", residues[i].SeqNumber, rmsdH.ToString("0.000"), rmsdE.ToString("0.000"));
                        // wRMSD.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", residues[i].SeqNumber, rmsdH.ToString("0.000"), rmsdH5.ToString("0.000"), rmsdH6.ToString("0.000"), rmsdE.ToString("0.000"));
                    }
                    wRMSD.Close();

                    times.Add(("Calculate RMSDs for whole chain", DateTime.Now.Subtract(stamp)));
                    stamp = DateTime.Now;
                }
                #endregion


                #region Calculating line segments corresponding to helices in template and query protein.

                List<double>[] dump;
                if (forceCalculateVectors || !tSSEs.All(sse => sse is SseInSpace))
                {
                    tSSEsInSpace = LibAnnotation.SSEsAsLineSegments_GeomVersion(tProtein.GetChain(templateChainID_), tSSEs, out dump);
                }
                else
                {
                    tSSEsInSpace = tSSEs.Select(sse => sse as SseInSpace).ToList();
                }

                if (forceCalculateVectors || !qSSEs.All(sse => sse is SseInSpace))
                {
                    qSSEsInSpace = LibAnnotation.SSEsAsLineSegments_GeomVersion(qProtein.GetChain(queryChainID), qSSEs, out dump);
                }
                else
                {
                    qSSEsInSpace = qSSEs.Select(sse => sse as SseInSpace).ToList();
                }

                times.Add(("Calculate line segments", DateTime.Now.Subtract(stamp)));
                stamp = DateTime.Now;
                #endregion


                #region Matching.

                List<Residue> tResiduesForAlignment = tProtein.GetChain(templateChainID_).GetResidues().ToList();
                List<Residue> qResiduesForAlignment = qProtein.GetChain(queryChainID).GetResidues().ToList();
                List<(int?, int?, double)> alignment = LibAnnotation.AlignResidues_DynProg(tResiduesForAlignment, qResiduesForAlignment);
                (List<int>, List<int>) pos = LibAnnotation.AlignmentToPositions(alignment);
                // double R = LibAnnotation.CharacteristicDistanceOfAlignment(tResiduesForAlignment, qResiduesForAlignment);
                // Console.WriteLine($"Characteristic distance R = {R}");

                Func<SseInSpace, SseInSpace, double> metricToMinimize;
                switch (metricMethod)
                {
                    case Setting.MetricMethod.No3:
                        metricToMinimize = LibAnnotation.MetricNo3Pos;
                        break;
                    case Setting.MetricMethod.No7:
                        metricToMinimize = (s, t) =>
                            LibAnnotation.MetricNo7Pos(s, t, alignment, LibAnnotation.DictResiToAli(tResiduesForAlignment, pos.Item1), LibAnnotation.DictResiToAli(qResiduesForAlignment, pos.Item2));
                        break;
                    case Setting.MetricMethod.No8:
                        Func<SseInSpace, SseInSpace, double> metric3 = LibAnnotation.MetricNo3Pos;
                        Func<SseInSpace, SseInSpace, double> metric7 = (s, t) => LibAnnotation.MetricNo7Pos(s, t, alignment, LibAnnotation.DictResiToAli(tResiduesForAlignment, pos.Item1), LibAnnotation.DictResiToAli(qResiduesForAlignment, pos.Item2));
                        double xx = 0.5;
                        metricToMinimize = (s, t) =>
                            (1 - xx) * metric3(s, t) + xx * metric7(s, t) + LibAnnotation.LengthDiffPenalty(s, t);
                        break;
                    default:
                        throw new Exception("Unknown MetricMethod: " + metricMethod);
                }


                Lib.WriteLineDebug("Template connections: {0}", tConnectivity.Select(t => tSSEs[t.Item1].Label + "-" + tSSEs[t.Item2].Label).EnumerateWithCommas());
                Lib.WriteLineDebug("Query connections: {0}", qConnectivity.Select(t => qSSEs[t.Item1].Label + "-" + qSSEs[t.Item2].Label).EnumerateWithCommas());

                AnnotationContext context = new AnnotationContext(metricToMinimize, annotationTypeFitting,
                    skipTemplatePenalty, skipCandidatePenalty,
                    tSSEsInSpace, qSSEsInSpace);
                context.InitializeTemplateConnectivity(tConnectivity);
                context.InitializeCandidateConnectivity(qConnectivity);

                if (selectionMethod == Setting.SelectionMethod.None)
                {
                    double[,] metricMatrix = new double[context.Templates.Count(), context.Candidates.Count()];
                    for (int i = 0; i < context.Templates.Count(); i++)
                    {
                        for (int j = 0; j < context.Candidates.Count(); j++)
                        {
                            SseInSpace t = context.Templates[i];
                            SseInSpace c = context.Candidates[j];
                            metricMatrix[i, j] = annotationTypeFitting(t.Type, c.Type) ? metricToMinimize(t, c) : -metricToMinimize(t, c);
                        }
                    }
                    AnnotationHelper.PrintCrossMatrix(context, metricMatrix, "metric_matrix.tsv");
                    Environment.Exit(0);
                }

                Func<AnnotationContext, IAnnotator> createAnnotator;
                switch (selectionMethod)
                {
                    case Setting.SelectionMethod.DynProg:
                        createAnnotator = DynProgAnnotator.New;
                        break;
                    case Setting.SelectionMethod.BB:
                        if (alternativeMatching)
                            context = context.WithAlternativeTemplates(templateSSA.MergeableSSEs);
                        if (softMatching)
                            context = context.Ordered().Softened_New(maxGapForSoftMatching);
                        createAnnotator = BranchAndBoundAnnotator.New;
                        break;
                    case Setting.SelectionMethod.MOM:
                        if (momTimeoutSeconds.HasValue)
                        {
                            createAnnotator = cont => new FallbackAnnotator(
                                new MOMAnnotator(cont, softMatching),
                                momTimeoutSeconds.Value,
                                "MOM annotator timeout, falling back to DP annotator",
                                new DynProgAnnotator(cont));
                        }
                        else
                        {
                            createAnnotator = cont => new MOMAnnotator(cont, softMatching);
                        }
                        break;
                    case Setting.SelectionMethod.Combined:
                        if (alternativeMatching)
                            context = context.WithAlternativeTemplates(templateSSA.MergeableSSEs);
                        if (softMatching)
                            context = context.Ordered().Softened_New(maxGapForSoftMatching);
                        createAnnotator = CombinedAnnotator.New;
                        break;
                    default:
                        throw new NotImplementedException("Not implemented for this selection method: " + selectionMethod.ToString());
                }
                NiceAnnotatorWrapper annotator =
                    correctionsFile == null ?
                    new NiceAnnotatorWrapper(context, createAnnotator)
                    : new NiceAnnotatorWrapperWithCorrections(context, createAnnotator, correctionsFile, queryID, qProtein);

                DateTime t_BB_0 = DateTime.Now;
                List<SseInSpace> annotQHelicesInSpace = annotator.GetAnnotatedCandidates().ToList();
                Lib.WriteLineDebug("Annotated: {0}", annotQHelicesInSpace.EnumerateWithSeparators("\n\t"));
                List<(int, int, int)> annotQConnectivity = annotator.GetAnnotatedConnectivity(qConnectivity);
                IEnumerable<double> metricList = annotator.GetMetricList();
                IEnumerable<double> suspiciousnessList = annotator.GetSuspiciousnessList();
                foreach (var sse in annotator.GetAnnotatedCandidates())
                {
                    if (sse.IsNotFound())
                    {
                        rmsdLists_AllChains.Add(new List<double>());
                    }
                    else
                    {
                        List<double>[] outRmsdLists;
                        LibAnnotation.SSEsAsLineSegments_GeomVersion(qProtein.GetChain(queryChainID), new Sse[] { sse }.ToList(), out outRmsdLists);
                        rmsdLists_AllChains.Add(outRmsdLists[0]);
                    }
                }
                List<String> qSseSequences = qProtein.GetChain(queryChainID)
                    .GetResidues(annotQHelicesInSpace.Select(s => (s.Start - extraResiduesInSequence, s.End + extraResiduesInSequence)))
                    .Select(s => String.Concat(s.Select(r => r.ShortName))).ToList();


                //var ssePairs = Enumerable.Zip (tSSEsInSpace,annotQHelicesInSpace,(t,q)=>new {T=t,Q=q});

                annotQHelicesInSpace_AllChains.AddRange(annotQHelicesInSpace);
                metricList_AllChains.AddRange(metricList);
                suspiciousnessList_AllChains.AddRange(suspiciousnessList);

                annotQSseSequences_AllChains.AddRange(qSseSequences);
                // metric3List_AllChains.AddRange (ssePairs.Select (x => x.Q.IsNotFound () ? Double.NaN : metric3(x.T,x.Q)));
                // metric7List_AllChains.AddRange (ssePairs.Select (x => x.Q.IsNotFound () ? Double.NaN : metric7(x.T,x.Q)));

                // metric3List_AllChains.AddRange (annotator.SelectFromAnnotated ((t,q) => q.IsNotFound () ? Double.NaN : metric3(t,q)));
                // metric7List_AllChains.AddRange (annotator.SelectFromAnnotated ((t,q) => q.IsNotFound () ? Double.NaN : metric7(t,q)));
                if (queryChainIDs.Length == 1)
                {
                    IEnumerable<(string, string)> labelPairs = annotator.GetMatching().Select(t =>
                      (annotator.Context.Templates[t.Item1].Label, annotator.Context.Candidates[t.Item2].Label));
                    if (Lib.DoWriteDebug)
                    {
                        using (StreamWriter w = new StreamWriter(Path.Combine(Setting.Directory, "matching-" + templateID + "-" + queryID + ".tsv")))
                        {
                            w.WriteLine("{0}\t{1}", templateID, queryID);
                            foreach (var t in labelPairs)
                            {
                                w.WriteLine("{0}\t{1}", t.Item1, t.Item2);
                            }
                        }
                    }
                }
                else
                {
                    Lib.WriteWarning("File matching-" + templateID + "-" + queryID + ".tsv was not produced, because annotating multiple chains.");
                }

                List<String> detQSseSequences = qProtein.GetChain(queryChainID)
                    .GetResidues(qSSEsInSpace.Select(s => (s.Start, s.End)))
                    .Select(s => String.Concat(s.Select(r => r.ShortName))).ToList();
                detQSsesInSpace_AllChains.AddRange(qSSEsInSpace);
                detQSseSequences_AllChains.AddRange(detQSseSequences);
                annotQConnectivity_AllChains.AddRange(annotQConnectivity);

                times.Add(("Matching", DateTime.Now - stamp));
                stamp = DateTime.Now;
                #endregion

            }
            // }
            #endregion


            if (Setting.FILTER_OUTPUT_BY_LABEL)
            {
                int[] indices = annotQHelicesInSpace_AllChains.IndicesWhere(sse => Setting.OUTPUT_ONLY_THESE_LABELS.Contains(sse.Label)).ToArray();
                annotQHelicesInSpace_AllChains = indices.Select(i => annotQHelicesInSpace_AllChains[i]).ToList();
                suspiciousnessList_AllChains = indices.Select(i => suspiciousnessList_AllChains[i]).ToList();
                rmsdLists_AllChains = indices.Select(i => rmsdLists_AllChains[i]).ToList();
                annotQSseSequences_AllChains = indices.Select(i => annotQSseSequences_AllChains[i]).ToList();
                metric3List_AllChains = indices.Select(i => metric3List_AllChains[i]).ToList();
                metric7List_AllChains = indices.Select(i => metric7List_AllChains[i]).ToList();
                annotQConnectivity_AllChains = new Lib.Shuffler(indices).UpdateIndices(annotQConnectivity_AllChains).ToList();
            }


            #region Writing values of metric.
            Lib.WriteInColor(ConsoleColor.Yellow, "Values of metric:\n");
            double totalMetric = metricList_AllChains.Where(x => !Double.IsNaN(x) && !Double.IsInfinity(x)).Sum();
            int totalFound = annotQHelicesInSpace_AllChains.Count(sse => !sse.IsNotFound());
            for (int i = 0; i < annotQHelicesInSpace_AllChains.Count; i++)
            {
                SseInSpace sse = annotQHelicesInSpace_AllChains[i];
                double metric = metricList_AllChains[i];
                Console.WriteLine("    {0,-5}: {1,8:0.00} {2}",
                    sse.Label,
                    sse.IsNotFound() ? Double.PositiveInfinity : metric,
                    sse.IsNotFound() ? "(not found)" : "");
            }
            Console.WriteLine("    Total: {0,8:0.00}", totalMetric);
            Console.WriteLine();
            Lib.WriteInColor(ConsoleColor.Yellow, "Annotated SSEs: {0} of {1}\n", totalFound, templateSSA.SSEs.Count);
            Console.WriteLine();
            #endregion


            #region Output of chosen SSEs into a file.
            String comment = "Automatic annotation for " + queryID + " based on " + templateID + " template.\nProgram was called with these parameters: "
                + String.Concat(args.Select(x => x + " ")) + "\nTotal value of used metric: " + totalMetric.ToString("0.00");
            if (Setting.JSON_OUTPUT)
            {
                var extras = Lib.DoWriteDebug ?
                    new Dictionary<string, IEnumerable<object>> {
                        {LibAnnotation.JsNames.METRIC, metricList_AllChains.Select (x => x as object)},
                        {LibAnnotation.JsNames.SUSPICIOUSNESS, suspiciousnessList_AllChains.Select (x => x as object)},
                        {LibAnnotation.JsNames.LIST_RMSD, rmsdLists_AllChains.Select (l => l as object)},
                        {LibAnnotation.JsNames.COUNT_RMSD, rmsdLists_AllChains.Select (l => l.Count() as object)},
                        {LibAnnotation.JsNames.FIRST_RMSD, rmsdLists_AllChains.Select (l => l.FirstOrDefault () as object)},
                        {LibAnnotation.JsNames.LAST_RMSD, rmsdLists_AllChains.Select (l => l.LastOrDefault () as object)},
                        {LibAnnotation.JsNames.AVG_RMSD, rmsdLists_AllChains.Select (l => l.Count()>=1 ? l.Average () : 0.0 as object)},
                        {LibAnnotation.JsNames.MAX_RMSD, rmsdLists_AllChains.Select (l => l.Count()>=1 ? l.Max () : 0.0 as object)},
                        {LibAnnotation.JsNames.ARG_MAX_RMSD, rmsdLists_AllChains.Select (l => l.Count()>=1 ? l.ArgMax () : -1 as object)},
                        {LibAnnotation.JsNames.MAX_INTERNAL_RMSD, rmsdLists_AllChains.Select (l => l.Count()>=3 ? l.Take (l.Count()-1).Skip (1).Max () : 0.0 as object)},
						// {"metric3_value", metric3List_AllChains.Select (x => x as object)},
						// {"metric7_value", metric7List_AllChains.Select (x => x as object)},
						{LibAnnotation.JsNames.SEQUENCE, annotQSseSequences_AllChains}
                    }
                    : new Dictionary<string, IEnumerable<object>> {
                        {LibAnnotation.JsNames.METRIC, metricList_AllChains.Select (x => x as object)},
                        {LibAnnotation.JsNames.SEQUENCE, annotQSseSequences_AllChains}
                    };
                LibAnnotation.WriteAnnotationFile_Json(fileQueryAnnotatedHelices, queryID,
                    annotQHelicesInSpace_AllChains,
                    extras,
                    annotQConnectivity_AllChains, // null, //TODO output beta_connectivity
                    querySSA.HBonds,
                    rotationMatrix,
                    comment);
            }
            else
            {
                LibAnnotation.WriteAnnotationFile(fileQueryAnnotatedHelices, annotQHelicesInSpace_AllChains, comment);
            }
            #endregion


            #region Running PyMOL to create .pse file.
            if (createPymolSession)
            {
                if (!Setting.JSON_OUTPUT)
                {
                    throw new NotImplementedException("Creating PyMOL session is implemented only for JSON output.");
                }
                Lib.WriteInColor(ConsoleColor.Yellow, "Running PyMOL:\n");
                if (!Lib.RunPyMOLScriptWithCommandLineArguments(config.PymolExecutable, config.PymolScriptSession, new string[]{
                    "cif",
                    Setting.Directory,
                    $"{templateID},{templateChainID_},{FormatRanges(templateDomainRanges)}",
                    $"{queryID},{queryChainIDs.Last()},{FormatRanges(queryDomainRanges)}",
                    }))
                    return -1;
            }

            times.Add(("Create PyMOL session", DateTime.Now.Subtract(stamp)));
            stamp = DateTime.Now;
            #endregion


            #region Write information about times.
            if (Lib.DoWriteDebug)
            {
                PrintTimes(times, t0);
            }
            #endregion


            return 0;
        }


        private static Protein ReadProteinFromFile(string filename, string chainId, IEnumerable<(int, int)> resSeqRanges)
        {
            try
            {
                Lib.WriteInColor(ConsoleColor.Yellow, "Loading structure:  {0}, chain {1}, residues {2}\n", filename, chainId, FormatRanges(resSeqRanges));
                Protein p;
                (int, int)[] resSeqRangesArray = resSeqRanges.Select(tup => (tup.Item1, tup.Item2)).ToArray();
                p = SecStrAnnotator2.CifWrapperForSecStrAnnot1.ProteinFromCifFile(filename, chainId, resSeqRangesArray);

                // Check entity type and emptiness of the structure
                if (p.GetChains().Any())
                {
                    string entityType = p.GetChains().First().GuessEntityType();
                    if (entityType != "protein")
                    {
                        Lib.WriteWarning($"Loaded structure is of type '{entityType}'");
                    }
                }
                else
                {
                    Lib.WriteWarning("Loaded structure contains no atoms.");
                }

                p = p.KeepOnlyNormalResidues(true);

                if (!p.GetChains().Any())
                {
                    Lib.WriteWarning("Loaded structure contains no normal residues.");
                    // Lib.WriteErrorAndExit("Loaded structure contains no normal residues.");
                }

                p = p.KeepOnlyOneAlternativeLocation();

                return p;
            }
            catch (IOException)
            {
                Lib.WriteErrorAndExit("Could not open \"" + filename + "\".");
                throw new Exception();
            }
        }


        private static void PrintTimes(List<(String, TimeSpan)> times, DateTime t0)
        {
            Console.WriteLine();
            Lib.WriteInColor(ConsoleColor.Yellow, "Times [miliseconds]:\n");
            times.ForEach(x => Console.WriteLine("> {0,6}   {1}", x.Item2.TotalMilliseconds.ToString("0"), x.Item1));
            Console.WriteLine("  ----------");
            Console.WriteLine("> {0}   Total", String.Format("{0,6}", DateTime.Now.Subtract(t0).TotalMilliseconds.ToString("0")));

        }


        private static List<(int, int)> ParseRanges(String rangeString)
        {
            List<(int, int)> result = new List<(int, int)>();
            foreach (String subrange in rangeString.Split(','))
            {
                string[] fromTo = subrange.Split(':').ToArray();
                if (fromTo.Length != 2)
                    throw new FormatException("Range must have format i:j, input was \"" + subrange + "\"");
                int fromI;
                int toI;
                try
                {
                    fromI = fromTo[0] != "" ? int.Parse(fromTo[0]) : int.MinValue;
                }
                catch (FormatException)
                {
                    throw new FormatException("Could not parse \"" + fromTo[0] + "\" as integer");
                }
                try
                {
                    toI = fromTo[1] != "" ? int.Parse(fromTo[1]) : int.MaxValue;
                }
                catch (FormatException)
                {
                    throw new FormatException("Could not parse \"" + fromTo[1] + "\" as integer");
                }
                //toI = fromTo[1] != "" ? int.Parse (fromTo[1]) : int.MaxValue;
                result.Add((fromI, toI));
            }
            return result;
        }


        private static String FormatRanges(IEnumerable<(int, int)> ranges)
        {
            return ranges.Select(
                range =>
                (range.Item1 == int.MinValue ? "" : range.Item1.ToString())
                + ":"
                + (range.Item2 == int.MaxValue ? "" : range.Item2.ToString())
            ).EnumerateWithSeparators(",");
        }


        private static (string pdb, string chain, List<(int, int)> ranges) ParseDomainSpecification(String domain, string defaultChain = null)
        {
            String[] parts = domain.Split(new char[] { ',' }, 3);

            string pdb = parts[0];

            string chain = (parts.Length >= 2) ? parts[1] : defaultChain;

            String rangeString = (parts.Length >= 3) ? parts[2] : Setting.DEFAULT_DOMAIN_RANGES_STRING;
            var ranges = ParseRanges(rangeString);

            return (pdb, chain, ranges);
        }
    }
}

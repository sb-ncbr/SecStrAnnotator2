using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;

using Cif.Components;
using protein.Sses;
using protein.Libraries;

namespace protein.SecStrAssigning
{
    public class DsspSecStrAssigner : ISecStrAssigner
    {
        public Protein Protein { get; private set; }
        public String DSSPExecutable { get; private set; }
        public String RenumberedPDBFile { get; private set; }  // CIF file with label_ values in auth_ fields
        public String DSSPFile { get; private set; }
        public string[] ChainIDs { get; private set; }
        public SseType[] AcceptedSSETypes { get; private set; }

        public DsspSecStrAssigner(Protein protein, String dsspExecutable, String renumberedPDBFile, String dsspFile, IEnumerable<string> chainIDs, SseType[] acceptedSSETypes)
        {
            this.Protein = protein;
            this.DSSPExecutable = dsspExecutable;
            this.RenumberedPDBFile = renumberedPDBFile;
            this.DSSPFile = dsspFile;
            this.ChainIDs = chainIDs.ToArray();
            this.AcceptedSSETypes = acceptedSSETypes;
        }

        public SecStrAssignment GetSecStrAssignment()
        {
            SecStrAssignment result;
            this.Protein.SaveCif(RenumberedPDBFile, fillAuthFieldsWithLabelValues: true, fabulateOccupancyAndBFactor: true);
            Lib.WriteInColor(ConsoleColor.Yellow, "Running DSSP:\n");
            if (!Lib.RunDSSP(DSSPExecutable, RenumberedPDBFile, DSSPFile))
            {
                throw new SecStrAssignmentException(MethodBase.GetCurrentMethod().DeclaringType.Name + " failed.");
            }
            try
            {
                // IEnumerable<SSE> sses = LibAnnotation.ReadSSEsFromDSSP (DSSPFile, AcceptedSSETypes)
                //     .Where (x => ChainIDs.Contains (x.ChainID))
                //     .Select ((sse, i) => sse.RelabeledCopy (sse.Type.AsString() + i.ToString ()));
                // result = new SecStrAssignment(sses); 
                return ReadSSEsFromDSSPFile();
                // TODO add support for beta-connectivity
            }
            catch (Exception e) { throw new SecStrAssignmentException(MethodBase.GetCurrentMethod().DeclaringType.Name + " failed.", e); }
            // File.Delete(RenumberedPDBFile);
            // File.Delete(DSSPFile);
            return result;
        }

        public String GetDescription()
        {
            return "DSSP method for " + Path.GetFileName(RenumberedPDBFile) + " (accepted types: " + Lib.EnumerateWithCommas(AcceptedSSETypes) + ")";
        }

        private SecStrAssignment ReadSSEsFromDSSPFile()
        {
            try
            {
                using (StreamReader r = new StreamReader(this.DSSPFile))
                {
                    SecStrAssignment assignment = ReadSSEsFromDSSPStreamReader(r);
                    return assignment;
                }
            }
            catch (FileNotFoundException e)
            {
                Lib.WriteError("Could not open \"{0}\".", this.DSSPFile);
                throw e;
            }
        }

        /** Reads DSSP output and returns a list of SSEs found in the protein structure.
		    Returns only SSEs of accepted types (G = 3_10 helix, H = alpha helix, I = pi helix...) */
        private SecStrAssignment ReadSSEsFromDSSPStreamReader(StreamReader reader)
        {
            List<Sse> SSEs = new List<Sse>();
            Dictionary<char, int?> dictSheetId = new Dictionary<char, int?> { { ' ', null } };
            var sse2ladders = new List<HashSet<char>>();
            int sheetIdCounter = 1;
            // TODO read beta-connections

            bool reading = false;
            String line;
            int counter = 0;
            char lastStrand1Char = ' ';
            char lastStrand2Char = ' ';
            char strand1Char = ' ';
            char strand2Char = ' ';
            while (!reader.EndOfStream)
            {
                line = reader.ReadLine();
                if (reading)
                {
                    String resSeqStr = line.Substring(5, 5);
                    // Lib.WriteWarning("resSeqStr: '{0}'", resSeqStr);
                    if (!resSeqStr.TrimStart().Equals(""))
                    {
                        int resSeq = Int32.Parse(resSeqStr);
                        // string chainID = line [11].ToString();
                        string chainID = line.Substring(153, 10).Trim();
                        // Lib.WriteWarning("chainID: '{0}'", chainID);
                        SseType type = LibSseTypes.Type(line[16].ToString());
                        lastStrand1Char = (strand1Char != ' ') ? strand1Char : lastStrand1Char;
                        lastStrand2Char = (strand2Char != ' ') ? strand2Char : lastStrand2Char;
                        strand1Char = line[23];
                        strand2Char = line[24];
                        char sheetChar = line[33];
                        int? sheetId = dictSheetId.GetOrAssignNext(sheetChar, ref sheetIdCounter);
                        if (this.AcceptedSSETypes.Contains(type))
                        {
                            Sse last = (SSEs.Count != 0) ? SSEs[SSEs.Count - 1] : null;
                            if (last != null && last.ChainID == chainID && last.End == resSeq - 1 && last.Type == type && sheetId == last.SheetId
                                    && DSSPStrandContinues(lastStrand1Char, lastStrand2Char, strand1Char, strand2Char))
                            {
                                last.End = resSeq;
                            }
                            else
                            {
                                Sse newSSE = new Sse(type.AsString() + counter.ToString(), chainID, resSeq, resSeq, type, sheetId);
                                SSEs.Add(newSSE);
                                sse2ladders.Add(new HashSet<char>());
                                last = newSSE;
                                Lib.WriteLineDebug("{0}, sheet char = {1}, sheet id = {2}", newSSE.Label, sheetChar, newSSE.SheetId);
                                counter++;
                            }
                            sse2ladders[counter - 1].Add(strand1Char);
                            sse2ladders[counter - 1].Add(strand2Char);
                        }
                    }
                }
                if (line.TrimStart()[0] == '#')
                    reading = true;
            }
            Dictionary<char, List<int>> ladder2sses = new Dictionary<char, List<int>>();
            for (int i = 0; i < counter; i++)
            {
                foreach (char strandChar in sse2ladders[i])
                {
                    if (strandChar != ' ')
                    {
                        if (!ladder2sses.ContainsKey(strandChar))
                        {
                            ladder2sses[strandChar] = new List<int>(2);
                        }
                        ladder2sses[strandChar].Add(i);
                    }
                }
            }
            List<(int, int, int)> betaConnectivity = new List<(int, int, int)>();
            foreach (char strandChar in ladder2sses.Keys)
            {
                var ssesHere = ladder2sses[strandChar];
                if (ssesHere.Count == 2)
                {
                    int orientation = Char.IsUpper(strandChar) ? -1 : Char.IsLower(strandChar) ? 1 : 0;
                    if (orientation == 0)
                    {
                        throw new Exception($"Invalid beta bridge label '{strandChar}', must be an uppercase or lowercase letter.");
                    }
                    int i = ssesHere.Min();
                    int j = ssesHere.Max();
                    betaConnectivity.Add((i, j, orientation));
                    Lib.WriteLineDebug($"beta  {SSEs[i].Label}  {SSEs[j].Label}  {orientation}");
                }
                else
                {
                    throw new Exception($"Ladder {strandChar} appears to connect {ssesHere.Count} strands (should always connect exactly 2 strands).");
                }
            }
            SecStrAssigning.SecStrAssignment assignment = new SecStrAssigning.SecStrAssignment(SSEs);
            assignment.Connectivity = betaConnectivity;
            Lib.WriteLineDebug("betaConnectivity: " + betaConnectivity.EnumerateWithCommas());
            return assignment;
        }

        private static bool DSSPStrandContinues(char lastStrand1, char lastStrand2, char currentStrand1, char currentStrand2)
        {
            return (lastStrand1 == currentStrand1 && currentStrand1 != ' ')
            || (lastStrand2 == currentStrand2 && currentStrand2 != ' ')
            || (lastStrand1 == currentStrand2 && currentStrand2 != ' ')
            || (lastStrand2 == currentStrand1 && currentStrand1 != ' ')
            || (currentStrand1 == ' ' && currentStrand2 == ' ')
            || (lastStrand1 == ' ' && lastStrand2 == ' ');
        }
    }
}


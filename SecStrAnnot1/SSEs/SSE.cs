using System;
using System.Collections.Generic;
using System.Linq;

using Cif.Components;
using protein.Libraries;

namespace protein.SSEs
{
    public class SSE : IComparable<SSE>
    {
        // An object of this class represents a secondary structure element (SSE) in a protein.

        public const int NOT_FOUND_START = 0;
        public const int NOT_FOUND_END = 0;
        public const string NOT_FOUND_CHAIN = "?";
        // public const SSEType NOT_FOUND_TYPE = SSEType.NOT_FOUND_TYPE;

        // public const char MIXED_HELIX_TYPE = 'h';
        // public const char HELIX_G_TYPE = 'G';
        // public const char HELIX_H_TYPE = 'H';
        // public const char HELIX_I_TYPE = 'I';
        // public const char SHEET_TYPE = 'E';
        // public const char ISOLATED_BETA_BRIDGE_TYPE = 'B';
        // public const char TURN_C7_TYPE = 'U';
        // public const char WIGGLE_C7_TYPE = 'W';
        // public const char BULGE_CLASSIC_SHORT_SIDE_TYPE = 'n';
        // public const char BULGE_CLASSIC_LONG_SIDE_TYPE = 'N';
        // public const char BULGE_WIDE_SHORT_SIDE_TYPE = 'm';
        // public const char BULGE_WIDE_LONG_SIDE_TYPE = 'M';
        // public const char BULGE_ANTIPARALLEL33_SHORT_SIDE_TYPE = 'u'; // in 2qad chain B ~ resi 15 // "short" side is the one donating protons
        // public const char BULGE_ANTIPARALLEL33_LONG_SIDE_TYPE = 'U'; // "long" side is the one accepting protons
        // public const char BULGE_ANTIPARALLEL22_SHORT_SIDE_TYPE = 't'; // in 1gei ~ resi 13 // "short" side is the one donating protons
        // public const char BULGE_ANTIPARALLEL22_LONG_SIDE_TYPE = 'T'; // "long" side is the one accepting protons
        // public const char BULGE_ANTIPARALLEL15_SHORT_SIDE_TYPE = 's'; // in 1gjm ~ resi 94
        // public const char BULGE_ANTIPARALLEL15_LONG_SIDE_TYPE = 'S'; 
        // public const char BULGE_ANTIPARALLEL23_SHORT_SIDE_TYPE = 'o'; // in 3dbg ~ resi 452
        // public const char BULGE_ANTIPARALLEL23_LONG_SIDE_TYPE = 'O'; 
        // public const char BULGE_PARALLEL14_SHORT_SIDE_TYPE = 'p';
        // public const char BULGE_PARALLEL14_LONG_SIDE_TYPE = 'P';
        // public const char BULGE_PARALLEL32_SHORT_SIDE_TYPE = 'q'; // in 3ruk ~ resi 38
        // public const char BULGE_PARALLEL32_LONG_SIDE_TYPE = 'Q';
        // public const char BULGE_PARALLEL13_SHORT_SIDE_TYPE = 'r'; // in 3dax ~ resi 35
        // public const char BULGE_PARALLEL13_LONG_SIDE_TYPE = 'R';
        // public const char BULGE_PARALLEL33_SHORT_SIDE_TYPE = 'l'; // in 3v8d ~ resi 69
        // public const char BULGE_PARALLEL33_LONG_SIDE_TYPE = 'L';
        // public static char[] ALL_HELIX_TYPES = { HELIX_H_TYPE, HELIX_G_TYPE, HELIX_I_TYPE, MIXED_HELIX_TYPE };
        // public static char[] ALL_SHEET_TYPES = { SHEET_TYPE, ISOLATED_BETA_BRIDGE_TYPE };

        public String Label { get; private set; } // Arbitrary name for the SSE.
        public string ChainID { get; private set; } // ID of the chain in which it is located.
        public int Start { get; set; } // Sequence number of the first residue.
        public int End { get; set; } // Sequence number of the last residue.
        public SSEType Type { get; set; } // Type of SSE according to DSSP abbreviations (H = alpha helix, G = 3_10 helix...)
        public int? SheetId { get; set; } // ID number of the beta-sheet that this beta-strand belong to. The value should be null for helices.
        public List<SSE> NestedSSEs { get; private set; }
        public String Comment { get; private set; } // Any string which is a comment for this SSE. 
        public String Color { get; set; } // Additional info about the color used for visualization (if assigned explicitly, else null).
        public string AuthChainID { get; private set; }
        public int? AuthStart { get; private set; }
        public string AuthStartInsCode { get; private set; }
        public int? AuthEnd { get; private set; }
        public string AuthEndInsCode { get; private set; }

        public SSE(String label, string chainID, int start, int end, SSEType type, int? sheetId)
        {
            Label = label;
            ChainID = chainID;
            Start = start;
            End = end;
            Type = type;
            SheetId = sheetId;
            NestedSSEs = null;
            Comment = null;
            Color = null;
        }
        public SSE(SSE orig)
        {
            Label = orig.Label;
            ChainID = orig.ChainID;
            Start = orig.Start;
            End = orig.End;
            Type = orig.Type;
            SheetId = orig.SheetId;
            NestedSSEs = orig.NestedSSEs;
            Comment = orig.Comment;
            Color = orig.Color;
            AuthChainID = orig.AuthChainID;
            AuthStart = orig.AuthStart;
            AuthStartInsCode = orig.AuthStartInsCode;
            AuthEnd = orig.AuthEnd;
            AuthEndInsCode = orig.AuthEndInsCode;
            // Console.WriteLine(this);
        }
        public static SSE NewNotFound(String label)
        {
            SSE result = new SSE(label, NOT_FOUND_CHAIN, NOT_FOUND_START, NOT_FOUND_END, SSEType.NOT_FOUND_TYPE, null);
            result.AddComment("Not found.");
            return result;
        }

        public bool IsNotFound()
        {
            return this.Type == SSEType.NOT_FOUND_TYPE;
        }

        public int Length()
        {
            return IsNotFound() ? 0 : End - Start + 1;
        }

        public override String ToString()
        {
            return "SSE " + Label + " in chain " + ChainID + " residues " + Start.ToString() + "-" + End.ToString() + " (type " + Type.ToString() + ")";
        }

        public override bool Equals(object obj)
        {
            if (obj is SSE)
            {
                SSE o = obj as SSE;
                return o.Label == this.Label && o.ChainID == this.ChainID && o.Start == this.Start && o.End == this.End && o.Type == this.Type;
            }
            else
            {
                return false;
            }
        }
        public override int GetHashCode()
        {
            return (Label?.GetHashCode() ?? 0) + ChainID.GetHashCode() + Start.GetHashCode() + End.GetHashCode() + Type.GetHashCode();
        }

        public static int Compare(SSE first, SSE second)
        {
            int result = first.ChainID.CompareTo(second.ChainID);
            if (result != 0)
                return result;
            result = first.Start.CompareTo(second.Start);
            if (result != 0)
                return result;
            return first.End.CompareTo(second.End);
        }
        public int CompareTo(SSE other)
        {
            return Compare(this, other);
        }

        /** Return an identical SSEInSpace, but with a different label. */
        public SSE RelabeledCopy(String newLabel)
        {
            SSE result = new SSE(this);
            result.Label = newLabel;
            return result;
        }

        /** Return an identical SSEInSpace, but with a different label. */
        public SSE RelabeledCopy(String newLabel, String newColor)
        {
            SSE result = new SSE(this);
            result.Label = newLabel;
            result.Color = newColor;
            return result;
        }

        public static SSE Join(params SSE[] sses)
        {
            if (sses.Length < 2)
            {
                throw new ArgumentException("The number of joined SSEs must be at least 2!");
            }
            var chainIDs = sses.Select(s => s.ChainID).Distinct().ToList();
            var sheetIDs = sses.Select(s => s.SheetId).Distinct().ToList();
            if (chainIDs.Count > 1)
            {
                throw new ArgumentException("Joined SSEs must be in the same chain!");
            }
            if (sheetIDs.Count > 1)
            {
                Lib.WriteWarning("Joining beta-strands with different sheet ID ({0})!", string.Join(", ", sheetIDs));
            }
            string newChainID = chainIDs[0];
            SSEType newType = sses.Select(sse => sse.Type).Aggregate<SSEType>((x, y) => Setting.JoiningTypeCombining(x, y) ?? SSEType.NOT_FOUND_TYPE);
            SSE first = sses[sses.Select(sse => sse.Start).ArgMin()];
            SSE last = sses[sses.Select(sse => sse.End).ArgMax()];
            SSE newSSE = new SSE(string.Join("+", sses.Select(s => s.Label ?? "")), newChainID, first.Start, last.End,
                newType, first.SheetId);
            newSSE.AuthChainID = first.AuthChainID;
            newSSE.AuthStart = first.AuthStart;
            newSSE.AuthStartInsCode = first.AuthStartInsCode;
            newSSE.AuthEnd = last.AuthEnd;
            newSSE.AuthEndInsCode = last.AuthEndInsCode;
            newSSE.AddComment("Created by joining " + sses.Count() + " SSEs: " + sses.Select(sse => sse.Label).EnumerateWithCommas() + ".");
            foreach (SSE sse in sses)
            {
                newSSE.AddNestedSSE(sse);
                if (sse.Comment != null)
                {
                    newSSE.AddComment(sse.Comment);
                }
            }
            return newSSE;
        }

        public void AddNestedSSE(SSE nested)
        {
            if (NestedSSEs == null)
                NestedSSEs = new List<SSE>();
            NestedSSEs.Add(nested);
        }

        /** Appends the specified string to this.Comment. */
        public void AddComment(String addedComment)
        {
            if (addedComment == null)
                return;
            else if (Comment == null)
                Comment = addedComment;
            else
                Comment = Comment + addedComment;
        }

        public bool IsSheet => Type.IsSheet();

        public bool IsHelix => Type.IsHelix();

        public void AddAuthFields(Protein protein)
        {
            // Console.WriteLine("AddAuthFields");
            try
            {
                Chain chain = protein.GetChain(this.ChainID);
                Residue res1 = chain.GetResidues().First(r => r.SeqNumber == this.Start);
                Residue res2 = chain.GetResidues().First(r => r.SeqNumber == this.End);
                this.AuthChainID = chain.AuthId;
                this.AuthStart = res1.AuthSeqNumber;
                this.AuthStartInsCode = res1.AuthInsertionCode;
                this.AuthEnd = res2.AuthSeqNumber;
                this.AuthEndInsCode = res2.AuthInsertionCode;
            }
            catch
            {
                Lib.WriteErrorAndExit($"Could not find starting or ending residue of SSE {ChainID} {Start}-{End}");
            }
        }
    }
}


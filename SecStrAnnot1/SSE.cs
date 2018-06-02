using System;
using System.Collections.Generic;
using System.Linq;

namespace protein
{
	public class SSE :IComparable<SSE>
	{
		// An object of this class represents a secondary structure element (SSE) in a protein.

		public const int NOT_FOUND_START = 0;
		public const int NOT_FOUND_END = 0;
		public const char NOT_FOUND_TYPE = 'X';
		public const char NOT_FOUND_CHAIN = 'X';
		public const char MIXED_HELIX_TYPE = 'h';
		public const char HELIX_G_TYPE = 'G';
		public const char HELIX_H_TYPE = 'H';
		public const char HELIX_I_TYPE = 'I';
		public const char SHEET_TYPE = 'E';
		public const char ISOLATED_BETA_BRIDGE_TYPE = 'B';
		public const char TURN_C7_TYPE = 'U';
		public const char WIGGLE_C7_TYPE = 'W';
		public const char BULGE_CLASSIC_SHORT_SIDE_TYPE = 'n';
		public const char BULGE_CLASSIC_LONG_SIDE_TYPE = 'N';
		public const char BULGE_WIDE_SHORT_SIDE_TYPE = 'm';
		public const char BULGE_WIDE_LONG_SIDE_TYPE = 'M';
		public const char BULGE_ANTIPARALLEL22_SHORT_SIDE_TYPE = 't'; // in 1gei ~ resi 13 // "short" side is the one donating protons
		public const char BULGE_ANTIPARALLEL22_LONG_SIDE_TYPE = 'T'; // "long" side is the one accepting protons
		public const char BULGE_ANTIPARALLEL15_SHORT_SIDE_TYPE = 's'; // in 1gjm ~ resi 94
		public const char BULGE_ANTIPARALLEL15_LONG_SIDE_TYPE = 'S'; 
		public const char BULGE_ANTIPARALLEL23_SHORT_SIDE_TYPE = 'o'; // in 3dbg ~ resi 452
		public const char BULGE_ANTIPARALLEL23_LONG_SIDE_TYPE = 'O'; 
		public const char BULGE_PARALLEL14_SHORT_SIDE_TYPE = 'p';
		public const char BULGE_PARALLEL14_LONG_SIDE_TYPE = 'P';
		public const char BULGE_PARALLEL32_SHORT_SIDE_TYPE = 'q'; // in 3ruk ~ resi 38
		public const char BULGE_PARALLEL32_LONG_SIDE_TYPE = 'Q';
		public const char BULGE_PARALLEL13_SHORT_SIDE_TYPE = 'r'; // in 3dax ~ resi 35
		public const char BULGE_PARALLEL13_LONG_SIDE_TYPE = 'R';
		public const char BULGE_PARALLEL33_SHORT_SIDE_TYPE = 'l'; // in 3v8d ~ resi 69
		public const char BULGE_PARALLEL33_LONG_SIDE_TYPE = 'L';
		public static char[] ALL_HELIX_TYPES = { HELIX_H_TYPE, HELIX_G_TYPE, HELIX_I_TYPE, MIXED_HELIX_TYPE };
		public static char[] ALL_SHEET_TYPES = { SHEET_TYPE, ISOLATED_BETA_BRIDGE_TYPE };

		public String Label { get; private set;} // Arbitrary name for the SSE.
		public char ChainID { get; private set;} // ID of the chain in which it is located.
		public int Start{ get; set;} // Sequence number of the first residue.
		public int End{ get; set;} // Sequence number of the last residue.
		public char Type{ get; set;} // Type of SSE according to DSSP abbreviations (H = alpha helix, G = 3_10 helix...)
		public int? SheetId{get;set;} // ID number of the beta-sheet that this beta-strand belong to. The value should be null for helices.
		public List<SSE> NestedSSEs{get; private set;}
		public String Comment { get; private set;} // Any string which is a comment for this SSE. 

		public SSE (String label, char chainID, int start, int end, char type, int? sheetId)
		{
			Label = label;
			ChainID = chainID;
			Start = start;
			End = end;
			Type = type;
			SheetId = sheetId;
			NestedSSEs = null;
			Comment = null;
		}
		public SSE (SSE orig){
			Label = orig.Label;
			ChainID = orig.ChainID;
			Start = orig.Start;
			End = orig.End;
			Type = orig.Type;
			SheetId = orig.SheetId;
			NestedSSEs = orig.NestedSSEs;
			Comment = orig.Comment;
		}
		public static SSE NewNotFound(String label){
			SSE result= new SSE(label,NOT_FOUND_CHAIN,NOT_FOUND_START, NOT_FOUND_END,NOT_FOUND_TYPE,null);
			result.AddComment ("Not found.");
			return result;
		}

		public bool IsNotFound(){
			return this.Type == NOT_FOUND_TYPE;
		}

		public int Length() {
			return IsNotFound () ? 0 : End - Start + 1;
		}

		public override String ToString(){
			return "SSE " + Label + " in chain " + ChainID + " residues " + Start.ToString () + "-" + End.ToString () + " (type " + Type.ToString () + ")"; 
		}

		public override bool Equals (object obj)
		{
			if (obj is SSE) {
				SSE o = obj as SSE;
				return o.Label == this.Label && o.ChainID == this.ChainID && o.Start == this.Start && o.End == this.End && o.Type == this.Type;
			} else {
				return false;
			}
		}
		public override int GetHashCode ()
		{
			return (Label?.GetHashCode ()??0) + ChainID.GetHashCode () + Start.GetHashCode () + End.GetHashCode () + Type.GetHashCode ();
		}

		public static int Compare(SSE first, SSE second){
			int result = first.ChainID.CompareTo (second.ChainID);
			if (result != 0)
				return result;
			result = first.Start.CompareTo (second.Start);
			if (result != 0)
				return result;
			return first.End.CompareTo (second.End);
		}
		public int CompareTo(SSE other){
			return Compare (this, other);
		}

		/** Return an identical SSEInSpace, but with a different label. */
		public SSE RelabeledCopy (String newLabel)
		{
			SSE result = new SSE (this);
			result.Label = newLabel;
			return result;
		}

		public static SSE Join(SSE first, SSE second, String comment){
			if (first.ChainID!=second.ChainID) throw new ArgumentException("Joined SSEs must be in the same chain!");
			if (first.Start>second.End) throw new ArgumentException("The first SSE must start before the second SSE's end!");
			if (first.SheetId != second.SheetId)
				Lib.WriteWarning ("Joining two beta-strands with different sheet ID ({0} and {1})!", first.SheetId, second.SheetId);
			SSE result = new SSE (first.Label + "_" + second.Label, first.ChainID, first.Start, second.End, first.Type == second.Type ? first.Type : SSE.MIXED_HELIX_TYPE,first.SheetId);
			result.AddComment (first.Comment);
			result.AddComment (second.Comment);
			result.AddComment ("Created by joining " + first.Start + "-" + first.End + " (type " + first.Type + ") and "
				+ second.Start + "-" + second.End + " (type " + second.Type + ")" + (comment == null ? "." : ", " + comment+"."));
			result.AddNestedSSE (first);
			result.AddNestedSSE (second);
			return result;
		}

		public void AddNestedSSE(SSE nested){
			if (NestedSSEs == null)
				NestedSSEs = new List<SSE> ();
			NestedSSEs.Add (nested);
		}

		/** Appends the specified string to this.Comment. */
		public void AddComment(String addedComment){
			if (addedComment == null)
				return;
			else if (Comment == null)
				Comment = addedComment;
			else
				Comment = Comment + addedComment;
		}

		public bool IsSheet{get{ return ALL_SHEET_TYPES.Contains(this.Type); }}
		public bool IsHelix{get{ return ALL_HELIX_TYPES.Contains(this.Type); }}


	}
}


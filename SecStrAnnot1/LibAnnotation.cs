using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Collections;
using Cif.Components;
using Cif.Tables;

namespace protein
{
	public static class LibAnnotation
	{
		public static bool alternativeLocationWarningPrinted = false;
		public const char SEPARATOR = '\t';
		public const char COMMENT_SIGN_READ = '#';
		public const String COMMENT_SIGN_WRITE = "# ";
		public class JsNames{
			public const string LABEL = "label";
			public const string CHAIN_ID = "chain_id";
			public const string START_RESI = "start";
			public const string END_RESI = "end";
			public const string TYPE = "type";
			public const string SHEET_ID = "sheet_id";
			public const string START_VECTOR = "start_vector";
			public const string END_VECTOR = "end_vector";
			public const string COMMENT = "comment";
			public const string COLOR = "color";
			public const string METRIC = "metric_value";
			public const string SUSPICIOUSNESS = "suspiciousness";

			public const string LIST_RMSD = "list_rmsd";
			public const string COUNT_RMSD = "count_rmsd";
			public const string FIRST_RMSD = "first_rmsd";
			public const string LAST_RMSD = "last_rmsd";
			public const string AVG_RMSD = "avg_rmsd";
			public const string MAX_RMSD = "max_rmsd";
			public const string ARG_MAX_RMSD = "arg_max_rmsd";
			public const string MAX_INTERNAL_RMSD = "max_internal_rmsd";

			public const string SEQUENCE = "sequence";
			public const string NESTED_SSES = "nested_sses";
			public const string COMMENT_PER_ENTRY = "comment";
			public const string FOUND_COUNT = "found";
			public const string TOTAL_METRIC = "total_metric_value";
			public const string SSES = "secondary_structure_elements";
			public const string BETA_CONNECTIVITY = "beta_connectivity";
			public const string H_BONDS = "hydrogen_bonds";
			public const string SSE_MERGING = "sse_merging";
			//public const string MERGED_LABEL = "label";
			//public const string FIRST_MERGED = "first";
			//public const string LAST_MERGED = "last";
		}

		const bool WRITE_EXTRAS = true;

		const bool WRITE_VECTORS = WRITE_EXTRAS;
		const bool WRITE_NESTED_SSES = WRITE_EXTRAS;
		const bool WRITE_METRIC = WRITE_EXTRAS;
		const bool WRITE_SUSPICIOUSNESS = WRITE_EXTRAS;
		const bool WRITE_SEQUENCE = WRITE_EXTRAS;
		const bool WRITE_COMMENTS_PER_ENTRY = WRITE_EXTRAS;

		const bool OUTPUT_NOT_FOUND_SSES = false;

		const int DEFAULT_DOUBLE_DIGITS = 6;

		const double SUSPICIOUSNESS_THRESHOLD = 1.0; // annotation of an SSE is suspicious, if annotated/(alternative+annotated)>= SUSPICIOUSNESS_THRESHOLD (annotated = really used metric, alternative = second best metric)



		/** Reads annotation data from file fileName in CSV format.
			CSV format: 1. column = Label, 2. column = ChainID, 3. column = Start, 4. column = End, 5. column = Type. 
						Columns are separated by SEPARATOR and comments are introduced by COMMENT_SIGN_READ. */
		public static List<SSE> ReadAnnotationFile(String fileName){
			List<List<String>> table;
			List<SSE> result = new List<SSE> ();
			HashSet<String> labels = new HashSet<String> (); // to check for duplicities
			using (StreamReader r = new StreamReader (fileName)) {
				table = Lib.ReadCSV (r, 5, SEPARATOR, COMMENT_SIGN_READ);
			}
			foreach (List<String> sse in table) {
				string chainID;
				int start;
				int end;
				char type;
				if (!labels.Add (sse [0]))
					Lib.WriteWarning ("Multiple definition of \"{0}\" in file \"{1}\".", sse [0], fileName);
				chainID = sse[1];
				try {
					start = Convert.ToInt32 (sse [2]);
				} catch (FormatException f) {
					throw new FormatException ("ReadAnnotationFile: Wrong format. Start must be an integer, but \"" + sse [2] + "\" was found in file \"\"+fileName+\"\".", f);
				}
				try {
					end = Convert.ToInt32 (sse [3]);
				} catch (FormatException f) {
					throw new FormatException ("ReadAnnotationFile: Wrong format. End must be an integer, but \"" + sse [2] + "\" was found in file \"" + fileName + "\".", f);
				}
				if (sse [4].Length == 1)
					type = sse [4] [0];
				else
					throw new FormatException ("ReadAnnotationFile: Wrong format. ChainID must be 1 character, but \"" + sse [1] + "\" was found in file \"\"+fileName+\"\".");

				result.Add (new SSE (sse [0], chainID, start, end, type,null));
			}
			return result;
		}

		/** Write annotation data tofile fileName in CSV format.
			CSV format: 1. column = Label, 2. column = ChainID, 3. column = Start, 4. column = End, 5. column = Type. 
						Columns are separated by SEPARATOR and comments are introduced by COMMENT_SIGN_READ. */
		public static void WriteAnnotationFile(String fileName, IEnumerable<SSE> sses, String comment){
			List<List<String>> table = sses.Where(x=>x!=null).Select (x => new List<String> {
				x.Label,
				x.ChainID.ToString (),
				x.Start.ToString (),
				x.End.ToString (),
				x.Type.ToString (),
				x.Comment
			}).ToList();
			using (StreamWriter w = new StreamWriter (fileName)) {
				Lib.WriteCSVWithComments (w, table, SEPARATOR, COMMENT_SIGN_WRITE, comment + "\n\nlabel\tchainID\tstart\tend\ttype");
			}
		}
		/** Write annotation data tofile fileName in CSV format.
			CSV format: 1. column = Label, 2. column = ChainID, 3. column = Start, 4. column = End, 5. column = Type, 6. column = Sequence. 
						Columns are separated by SEPARATOR and comments are introduced by COMMENT_SIGN_READ. */
		public static void WriteAnnotationFileWithSequences(String fileName, IEnumerable<Tuple<SSE,String>> sses, String comment){
			List<List<String>> table = sses.Where(t=>t!=null).Select (t => new List<String> {
				t.Item1.Label,
				t.Item1.ChainID.ToString (),
				t.Item1.Start.ToString (),
				t.Item1.End.ToString (),
				t.Item1.Type.ToString (),
				t.Item2,
				t.Item1.Comment
			}).ToList();
			using (StreamWriter w = new StreamWriter (fileName)) {
				Lib.WriteCSVWithComments (w, table, SEPARATOR, COMMENT_SIGN_WRITE, comment + "\n\nlabel\tchainID\tstart\tend\ttype\tsequence");
			}
		}

		public static JsonValue SSEToJson(SSE sse,IDictionary<string,IEnumerator<object>> extraEnumerators/*IEnumerator<double> metricEnumerator,IEnumerator<double> suspiciousEnumerator,IEnumerator<String> sequenceEnumerator*/){
			JsonValue elem = JsonValue.MakeObject ();
			elem [JsNames.LABEL] = new JsonValue (sse.Label);
			elem [JsNames.CHAIN_ID] = new JsonValue (sse.ChainID.ToString ());
			elem [JsNames.START_RESI] = new JsonValue (sse.Start);
			elem [JsNames.END_RESI] = new JsonValue (sse.End);
			elem [JsNames.TYPE] = new JsonValue (sse.Type.ToString ());
			if (sse.SheetId != null)
				elem [JsNames.SHEET_ID] = new JsonValue ((int)sse.SheetId);
			if (sse.Comment != null)
				elem [JsNames.COMMENT] = new JsonValue (sse.Comment);
			if (sse.Color != null)
				elem [JsNames.COLOR] = new JsonValue (sse.Color);
			if (WRITE_VECTORS && sse is SSEInSpace) {
				elem [JsNames.START_VECTOR] = new JsonValue ((sse as SSEInSpace).StartVector.Round (DEFAULT_DOUBLE_DIGITS).AsList ());
				elem [JsNames.END_VECTOR] = new JsonValue ((sse as SSEInSpace).EndVector.Round (DEFAULT_DOUBLE_DIGITS).AsList ());
			}
			if (extraEnumerators != null) {
				foreach (var kv in extraEnumerators) {
					if (!kv.Value.MoveNext ())
						throw new Exception ();
					if (kv.Value.Current is double)
						elem [kv.Key] = new JsonValue (RoundDouble((double)kv.Value.Current));
					else if (kv.Value.Current is int)
						elem [kv.Key] = new JsonValue (RoundDouble(1.0*(int)kv.Value.Current));
					else if (kv.Value.Current is string)
						elem [kv.Key] = new JsonValue ((string)kv.Value.Current);
					else if (kv.Value.Current is IEnumerable<double>)
						elem [kv.Key] = new JsonValue ((kv.Value.Current as IEnumerable<double>).Select (x=>x as object).ToList());
					else if (kv.Value.Current == null)
						elem [kv.Key] = new JsonValue ();
					else
						throw new Exception ("Unsupported type: " + kv.Value.Current.GetType ().ToString ());
				}
			}
			if (WRITE_NESTED_SSES && sse.NestedSSEs != null) {
				JsonValue nested = JsonValue.MakeList ();
				bool anyNested = false;
				foreach (SSE nestedSSE in sse.NestedSSEs) {
					if (Lib.DoWriteDebug || nestedSSE.Label != null) {
						nested.Add (SSEToJson (nestedSSE, null));
						anyNested = true;
					}
				}
				if (anyNested) {
					elem [JsNames.NESTED_SSES] = nested;
				}
			}
			return elem;
		}

		public static void WriteAnnotationFile_Json(String fileName, String name, 
			IEnumerable<SSE> sses, IDictionary<String,IEnumerable<object>> extras, List<Tuple<int,int,int>> betaConnectivity, List<Tuple<Residue,Residue>> hBonds, String comment){

			JsonValue ssesJson = JsonValue.MakeList ();
			Dictionary<String,IEnumerator<object>> extraEnumerators = extras?.ToDictionary (kv => kv.Key, kv => kv.Value.GetEnumerator ());

			foreach (SSE sse in sses.Where (x => x != null)) {
				JsonValue elem = SSEToJson (sse, extraEnumerators);
				if (OUTPUT_NOT_FOUND_SSES || !sse.IsNotFound ())
					ssesJson.Add (elem);
			}

			JsonValue json = JsonValue.MakeObject ();
			json [name] = JsonValue.MakeObject ();
			if (WRITE_COMMENTS_PER_ENTRY && comment != null)
				json [name] [JsNames.COMMENT_PER_ENTRY] = new JsonValue (comment?.Replace ('\n', ' ')?.Replace ('\t', ' '));
			json [name] [JsNames.FOUND_COUNT] = new JsonValue (sses.Count (sse => !sse.IsNotFound ()));
			if (WRITE_METRIC && extras!=null && extras.ContainsKey (JsNames.METRIC))
				json [name] [JsNames.TOTAL_METRIC] = new JsonValue (extras[JsNames.METRIC].Select (x=> RoundDouble((double)x)).Where (x=>!Double.IsInfinity (x)&& !Double.IsNaN (x)).Sum ());
			json [name] [JsNames.SSES] = ssesJson;
			if (Lib.DoWriteDebug && hBonds != null)
				json [name] [JsNames.H_BONDS] = HBondsToJson (hBonds);
			if (betaConnectivity != null)
				json [name] [JsNames.BETA_CONNECTIVITY] = BetaConnectivityToJson (betaConnectivity, sses);

			using (StreamWriter w = new StreamWriter (fileName)) {
				w.Write (json.ToString(MainClass.JSON_OUTPUT_MAX_INDENT_LEVEL));
			}
		}

		public static JsonValue HBondsToJson(List<Tuple<Residue,Residue>> hBonds){
			JsonValue hBondsJson = JsonValue.MakeList ();
			foreach (var bond in hBonds) {
				JsonValue bondJson = JsonValue.MakeList ();
				bondJson.Add (new JsonValue(bond.Item1.ChainId.ToString ()));
				bondJson.Add (new JsonValue(bond.Item1.SeqNumber));
				bondJson.Add (new JsonValue(bond.Item2.ChainId.ToString ()));
				bondJson.Add (new JsonValue(bond.Item2.SeqNumber));
				hBondsJson.Add (bondJson);
			}
			return hBondsJson;
		}

		public static JsonValue BetaConnectivityToJson(List<Tuple<int,int,int>> betaConnectivity, IEnumerable<SSE> sses){
			if (betaConnectivity == null)
				return new JsonValue ();
			List<String> labels = sses.Select (sse => sse.Label).ToList ();
			JsonValue connsJson = JsonValue.MakeList ();
			foreach (var conn in betaConnectivity) {
				JsonValue connJson = JsonValue.MakeList ();
				try{
					//connJson.Add (new JsonValue(conn.Item1));
					//connJson.Add (new JsonValue(conn.Item2));
					connJson.Add (new JsonValue(labels[conn.Item1]));
					connJson.Add (new JsonValue(labels[conn.Item2]));
					connJson.Add (new JsonValue(conn.Item3));
					connsJson.Add (connJson);
				} catch (ArgumentOutOfRangeException e){
					// do not output this conn
				}
			}
			return connsJson;
		}

		/* Takes the first entity in the file. */
		public static List<SSE> ReadAnnotationFile_Json(String fileName){
			String dump;
			List<Tuple<int,int,int>> dump2;
			List<Tuple<String,int,int>> dump3;
			return ReadAnnotationFile_Json (fileName, null, out dump, out dump2, out dump3, false);
		}

		/* If name==null, takes the first entity in the file. */
		public static List<SSE> ReadAnnotationFile_Json(String fileName, String name, out String comment, out List<Tuple<int,int,int>> betaConnectivity, out List<Tuple<String,int,int>> merging, bool takeTheOnlyEntryEvenIfNameDoesntMatch){
			List<SSE> result = new List<SSE> ();
			String str;
			using (StreamReader r = new StreamReader (fileName)) {
				str = r.ReadToEnd ();
			}
			JsonValue json;
			try {
				json = JsonValue.FromString (str);
			}catch (Exception){//TODO put actual exception type here
				throw new FormatException (fileName + " is not a valid JSON.");
			}
				
			if (json.Type != JsonType.Object)
				throw new FormatException (fileName+" is not a JSON object.");
			if (json.Count==0)
				throw new FormatException (fileName + " is an empty JSON object.");
			if (name == null)
				name = json.Key (0);
			if (!json.Contains (name)) {
				if (json.Count == 1) {
					Lib.WriteWarning(fileName + " does not contain entry '" + name + "', taking '" + json.Key(0) + "' instead.");
					name = json.Key(0);
				} else {
					throw new FormatException (fileName + " does not contain key " + name + ".");
				}
			}
			JsonValue entry = json [name];
			if (entry.Type!=JsonType.Object || !entry.Contains (JsNames.SSES))
				throw new FormatException (JsonLocationString (fileName,name)+" does not contain key \""+JsNames.SSES+"\".");
			if (entry.Contains (JsNames.COMMENT_PER_ENTRY))
				comment = entry [JsNames.COMMENT_PER_ENTRY].String;
			else
				comment = null;
			JsonValue sses = entry [JsNames.SSES];
			if (sses.Type != JsonType.List)
				throw new FormatException (JsonLocationString (fileName, name, JsNames.SSES) + " is not a JSON list.");
			for (int i = 0; i < sses.Count; i++) {
				JsonValue sse = sses [i];
				if (sse.Type != JsonType.Object)
					throw new FormatException (JsonLocationString (fileName,name,JsNames.SSES,i) + "  is not a JSON object.");
				try {
					String chainId = sse [JsNames.CHAIN_ID].String;
					String type = sse [JsNames.TYPE].String;
					if (type.Length != 1)
						throw new FormatException (JsonLocationString (fileName,name,JsNames.SSES,i,JsNames.TYPE) + " is not a single character.");
					int? sheetId=sse.Contains (JsNames.SHEET_ID) ? sse [JsNames.SHEET_ID].Int as int? : null;
					SSE s = new SSE (sse [JsNames.LABEL].String, chainId, sse [JsNames.START_RESI].Int, sse [JsNames.END_RESI].Int, type [0],sheetId);

					if (sse.Contains (JsNames.START_VECTOR) && sse.Contains (JsNames.END_VECTOR)) {
						JsonValue vec = sse [JsNames.START_VECTOR];
						Vector startVector = new Vector (vec [0].Double, vec [1].Double, vec [2].Double);
						vec = sse [JsNames.END_VECTOR];
						Vector endVector = new Vector (vec [0].Double, vec [1].Double, vec [2].Double);
						s = new SSEInSpace (s, startVector, endVector);
					}
					if (sse.Contains (JsNames.COMMENT))
						s.AddComment (sse [JsNames.COMMENT].String);
					if (sse.Contains (JsNames.COLOR))
						s.Color = sse[JsNames.COLOR].String;
					result.Add (s);
				} catch (JsonKeyNotFoundException e) {
					throw new FormatException ("Error in parsing JSON object " + JsonLocationString (fileName,name,JsNames.SSES,i) + ". " + "Does not contain key \"" + e.MissingKey + "\".");
				} catch (Exception e) {
					throw new FormatException ("Error in parsing JSON object " + JsonLocationString (fileName,name,JsNames.SSES,i) + ".", e);
				}
			}

			// Reading beta-connectivity
			bool labelDuplicates = result.Select (sse => sse.Label).Distinct ().Count () != result.Count;
			betaConnectivity = new List<Tuple<int, int,int>> ();
			if (entry.Contains (JsNames.BETA_CONNECTIVITY)) {
				JsonValue connections = entry [JsNames.BETA_CONNECTIVITY];
				if (connections.Type!=JsonType.List)
					throw new FormatException (JsonLocationString (fileName,name,JsNames.BETA_CONNECTIVITY)+" is not a JSON list.");
				JsonValue[] connsList = connections.Where (c => c.Type != JsonType.Null).ToArray (); //This null think must be here cuz the Json parser sucks.
				for (int i = 0; i < connsList.Length; i++) {
				//foreach (JsonValue connection in connections.Where (c=>c.Type!=JsonType.Null)/*This null think must be here cuz the Json parser sucks.*/) {
					JsonValue connection=connsList[i];
					if (connection.Type!=JsonType.List||connection.Count!=3)
						throw new FormatException (JsonLocationString (fileName,name,JsNames.BETA_CONNECTIVITY,i)+" is not a JSON list with three 3 elements.");
					if (connection.All (v => v.Type == JsonType.Int)) {
						betaConnectivity.Add (new Tuple<int, int, int> (connection [0].Int, connection [1].Int, connection [2].Int));
					} else if (connection [0].Type == JsonType.String && connection [1].Type == JsonType.String && connection [2].Type == JsonType.Int) {
						if (labelDuplicates)
							throw new FormatException (JsonLocationString (fileName,name,JsNames.SSES) + " contains label duplicates. Cannot unambigously read " + JsonLocationString (fileName,name,JsNames.BETA_CONNECTIVITY));
						int first = result.FindIndex (sse => sse.Label == connection [0].String);
						int last = result.FindIndex (sse => sse.Label == connection [1].String);
						if (first == -1)
							throw new FormatException (JsonLocationString (fileName,name,JsNames.BETA_CONNECTIVITY,i,0) + " is not a label found in " + JsNames.SSES + ".");
						if (last == -1)
							throw new FormatException (JsonLocationString (fileName,name,JsNames.BETA_CONNECTIVITY,i,1) + " is not a label found in " + JsNames.SSES + ".");
						betaConnectivity.Add (new Tuple<int, int, int> (first, last, connection [2].Int));
					} else {
						throw new FormatException (JsonLocationString (fileName,name,JsNames.BETA_CONNECTIVITY,i)+" must be a JSON list with three 3 elements (labels/indices of 2 SSEs and an integer for direction of the connection (1=parallel,-1=antiparallel)).");
					}
				}
			}

			// Reading mergeable SSEs
			merging = new List<Tuple<string, int, int>>();
			if (entry.Contains (JsNames.SSE_MERGING)) {
				JsonValue mergs = entry [JsNames.SSE_MERGING];
				if (mergs.Type!=JsonType.List)
					throw new FormatException (JsonLocationString (fileName,name,JsNames.SSE_MERGING)+" is not a JSON list.");
				JsonValue[] mergsList = mergs.Where (m => m.Type != JsonType.Null).ToArray (); //This null think must be here cuz the Json parser sucks.
				for (int i = 0; i < mergsList.Length; i++) {
					JsonValue merg = mergsList [i];
					if (merg.Type!=JsonType.List||merg.Count!=3)
						throw new FormatException (JsonLocationString (fileName,name,JsNames.SSE_MERGING,null)+" is not a JSON list with three 3 elements.");
					if (merg [0].Type == JsonType.String && merg [1].Type == JsonType.Int && merg [2].Type == JsonType.Int) {
						merging.Add (new Tuple<String, int, int> (merg [0].String, merg [1].Int, merg [2].Int));
					} else if (merg.All (m => m.Type == JsonType.String)) {
						if (labelDuplicates)
							throw new FormatException (JsonLocationString (fileName,name,JsNames.SSES) + " contains label duplicates. Cannot unambigously read " + JsonLocationString (fileName,name,JsNames.SSE_MERGING));
						int first = result.FindIndex (sse => sse.Label == merg [1].String);
						int last = result.FindIndex (sse => sse.Label == merg [2].String);
						if (first == -1)
							throw new FormatException (JsonLocationString (fileName,name,JsNames.SSE_MERGING,i,1) +" is not a label found in " + JsonLocationString (fileName,name,JsNames.SSES) + ".");
						if (last == -1)
							throw new FormatException (JsonLocationString (fileName,name,JsNames.SSE_MERGING,i,2) +" is not a label found in " + JsonLocationString (fileName,name,JsNames.SSES) + ".");
						merging.Add (new Tuple<String, int, int> (merg[0].String, first, last));
					} else {
						throw new FormatException (JsonLocationString (fileName,name,JsNames.SSE_MERGING,i)+" must be a JSON list with three 3 elements (1 label of merged SSE and 2 labels/indices of first and last SSE to be merge).");
					}
				}
			}
			return result;
		}

		/* null in keys will be formatted as unknown index _ */
		private static String JsonLocationString(string filename, params object[] keys){
			return filename + keys.Select (k => "[" + (k == null ? "_" : k is string ? "\"" + k + "\"" : k.ToString ()) + "]").EnumerateWithSeparators ("");
		}

		/*private static String JsonLocationValueString(object value, string filename, params object[] keys){
			return JsonLocationString(filename,keys) + ": " + (value == null ? "null" : value is string ? "\"" + value + "\"" : value.ToString ());
		}*/

		private static double RoundDouble(double x){
			return Math.Round (x, DEFAULT_DOUBLE_DIGITS);
		}

		/** Reads DSSP output and returns a list of SSEs found in the protein structure. */
		public static List<SSE> ReadSSEsFromDSSP(String fileName, char[] acceptedTypes){
			try {
				StreamReader r = new StreamReader (fileName);
				List<SSE> helices = ReadSSEsFromDSSP (r,acceptedTypes);
				r.Close ();
				return helices;
			} catch (FileNotFoundException e) {
				Lib.WriteError ("Could not open \"{0}\".", fileName);
				throw e;
			}
		}

		/** Reads DSSP output and returns a list of SSEs found in the protein structure.
			Returns only SSEs of accepted types (G = 3_10 helix, H = alpha helix, I = pi helix...) */
		public static List<SSE> ReadSSEsFromDSSP(StreamReader reader, char[] acceptedTypes){
			List<SSE> SSEs = new List<SSE> ();
			Dictionary<char,int?> dictSheetId = new Dictionary<char, int?> { { ' ',null } };
			int sheetIdCounter = 1;

			bool reading = false;
			String line;
			int counter = 0;
			char lastStrand1Char = ' ';
			char lastStrand2Char = ' ';
			char strand1Char = ' ';
			char strand2Char = ' ';
			while (!reader.EndOfStream) {
				line = reader.ReadLine ();
				if (reading) {
					String resSeqStr = line.Substring (6, 4);
					if (!resSeqStr.TrimStart ().Equals ("")) {
						int resSeq = Int32.Parse (resSeqStr);
						string chainID = line [11].ToString();
						char type = line [16];
						lastStrand1Char = (strand1Char != ' ') ? strand1Char : lastStrand1Char;
						lastStrand2Char = (strand2Char != ' ') ? strand2Char : lastStrand2Char;
						strand1Char = line [23];
						strand2Char = line [24];
						char sheetChar = line [33];
						int? sheetId = dictSheetId.GetOrAssignNext (sheetChar, ref sheetIdCounter);
						if (acceptedTypes.Contains(type)) {
							SSE last = (SSEs.Count != 0) ? SSEs[SSEs.Count-1] : null;
							if (last != null && last.ChainID == chainID && last.End == resSeq - 1 && last.Type == type && sheetId==last.SheetId && DSSPStrandContinues (lastStrand1Char,lastStrand2Char,strand1Char,strand2Char)) {
								last.End = resSeq;
							} else {
								counter++;
								SSE newSSE = new SSE (type.ToString () + counter.ToString (), chainID, resSeq, resSeq, type,sheetId);
								SSEs.Add (newSSE);
								//Lib.WriteLineDebug ("{0}, sheet char = {1}, sheet id = {2}", newSSE.Label,sheetChar,newSSE.SheetId);
							}
						}
					}
				}
				if (line.TrimStart () [0] == '#')
					reading = true;
			}
			return SSEs;
			//TODO Add checking the format of the file.
		}

		private static bool DSSPStrandContinues(char lastStrand1,char lastStrand2,char currentStrand1,char currentStrand2){
			return (lastStrand1 == currentStrand1 && currentStrand1 != ' ')
			|| (lastStrand2 == currentStrand2 && currentStrand2 != ' ')
			|| (lastStrand1 == currentStrand2 && currentStrand2 != ' ')
			|| (lastStrand2 == currentStrand1 && currentStrand1 != ' ')
			|| (currentStrand1 == ' ' && currentStrand2 == ' ')
			|| (lastStrand1 == ' ' && lastStrand2 == ' ');
		}

		public class JoiningCriteriaValues{
			public int Gap{ get; private set; }
			public double Angle{ get; private set; }
			public double Stagger{ get; private set; }
			public JoiningCriteriaValues(int gap, double angle,double stagger){
				Gap=gap;
				Angle=angle;
				Stagger=stagger;
			}
			public override String ToString(){
				return this.GetType ().Name+" [gap: " + Gap + ", angle: " + Angle.ToString ("0.00") + " deg, stagger: " + Stagger.ToString ("0.00") + " A]";
			}
		}

		/** A class encapsulating parameters for joining neighboring SSEs. */
		public class JoiningParameters{
			public char[] AcceptedTypes{ get; private set;}
			public double GapPenalty { get; private set;}
			public double AnglePenalty { get; private set;} // angles in degrees! 
			public double StaggerPenalty { get; private set; } // stagger is in Angstroms
			public int GapThreshold { get; private set;}
			public double PenaltyThreshold { get; private set;}
			public JoiningParameters(char[] acceptedTypes, double gapPenalty,double anglePenaltyPerDegree, double staggerPenalty, int gapThreshold,double penaltyThreshold){
				AcceptedTypes=acceptedTypes;
				GapPenalty=gapPenalty;
				AnglePenalty=anglePenaltyPerDegree;
				StaggerPenalty=staggerPenalty;
				GapThreshold=gapThreshold;
				PenaltyThreshold=penaltyThreshold;
			}
			public override String ToString (){
				return "AcceptedTypes = " + AcceptedTypes.Aggregate("", (x,y)=>x+y)
				+ ", GapPenalty = " + GapPenalty
				+ ", AnglePenalty = " + AnglePenalty
				+ ", GapThreshold = " + GapThreshold
				+ ", PenaltyThreshold = " + PenaltyThreshold;
			}
		}

		/** Joins neighboring SSEs if they fulfil some conditions. */
		public static List<SSEInSpace> JoinSSEs(List<SSEInSpace> sses, JoiningParameters parameters){
			List<SSEInSpace> result = new List<SSEInSpace> ();
			if (sses.Count == 0)
				return result;

			JoiningCriteriaValues criteria;

			SSEInSpace a = sses [0];
			foreach (SSEInSpace b in sses.Skip(1)) {
				if (ShouldJoin (a, b, parameters, out criteria)) {
					Console.WriteLine ("Joining {0} ({1}-{2}) and {3} ({4}-{5}) with gap {6} and angle {7}°.",a.Label,a.Start,a.End,b.Label,b.Start,b.End,b.Start - a.End - 1,180/Math.PI*Geom.AngleInRadians (a.EndVector - a.StartVector, b.EndVector - b.StartVector));
					a = SSEInSpace.Join (a, b, criteria.ToString ());
				} else {
					result.Add (a);
					a = b;
				}
			}
			result.Add (a);
			return result;
		}

		/** Says whether too SSEs should be joined. */
		private static bool ShouldJoin(SSEInSpace sse1, SSEInSpace sse2, JoiningParameters parameters, out JoiningCriteriaValues criteria){
			if (!parameters.AcceptedTypes.Contains (sse1.Type) || !parameters.AcceptedTypes.Contains (sse2.Type)) {
				criteria = null;
				return false;
			}
			if (sse1.ChainID != sse2.ChainID) {
				criteria = null;
				return false;
			}

			int gap = sse2.Start - sse1.End - 1;
			double angle = Geom.AngleInDegrees (sse1.LineSegment.Direction, sse2.LineSegment.Direction);
			double stagger = Stagger (sse1, sse2);
			criteria = new JoiningCriteriaValues (gap, angle, stagger);

			if (gap<=parameters.GapThreshold && gap * parameters.GapPenalty + angle * parameters.AnglePenalty + stagger * parameters.StaggerPenalty <= parameters.PenaltyThreshold)
				return true;
			else
				return false;
		}

		/** Joins neighboring SSEs if they fulfil some conditions. */
		public static List<SSE> JoinSSEs_GeomVersion(Protein p, List<SSE> sses, double rmsdLimit, Func<char,char,char?> typeCombining){
			List<SSE> result = new List<SSE> ();

			foreach (Chain chain in p.GetChains ()) {		
				List<SSE> chainSSEs = sses.Where (x => x.ChainID == chain.Id).OrderBy (x => x.Start).ToList ();

				List<Tuple<int,int>> ranges = new List<Tuple<int, int>> ();
				List<char?> checkTypes = new List<char?> ();
				for (int i = 0; i < chainSSEs.Count - 1; i++) {
					ranges.Add (new Tuple<int,int> (chainSSEs [i].End - 3, chainSSEs [i + 1].Start + 3));
					checkTypes.Add (typeCombining (chainSSEs [i].Type, chainSSEs [i + 1].Type));
				}

				double[] rmsds = new double[ranges.Count];
				List<GeometryCheckResult> checkResults = chain.GetResidues (ranges).Select ((l,i) => checkTypes[i]==null ? GeometryCheckResult.NOK: CheckGeometry (l, (char) checkTypes[i], rmsdLimit, out rmsds[i])).ToList ();

				if (chainSSEs.Count > 0) {
					result.Add (chainSSEs [0]);
				}
				for (int i = 1; i < chainSSEs.Count; i++) {
					if (checkResults [i - 1] == GeometryCheckResult.OK) {
						Console.WriteLine ("Joining {0} ({1}-{2}) and {3} ({4}-{5}).",
							result [result.Count - 1].Label, result [result.Count - 1].Start, result [result.Count - 1].End, chainSSEs [i].Label, chainSSEs [i].Start, chainSSEs [i].End);
						result [result.Count - 1] = SSE.Join (result [result.Count - 1], chainSSEs [i], "Max RMSD: "+rmsds[i-1].ToString ("0.000")+" A");
					} else {
						result.Add (chainSSEs [i]);
					}
				}
			}
			return result;
		}

		/** Calculates the "stagger" metric for two helices, which is a measure of their continuity. The value is in Angstroms. */
		private static double Stagger (SSEInSpace sse1, SSEInSpace sse2){
			Geom.Point A = new Geom.Point (sse1.StartVector);
			Geom.Point B = new Geom.Point (sse1.EndVector);
			Geom.Point C = new Geom.Point (sse2.StartVector);
			Geom.Point D = new Geom.Point (sse2.EndVector);
			Geom.Plane rho = new Geom.Plane (Geom.Middle (B,C), D.Vector-A.Vector);
			Geom.Point P = Geom.Intersection (new Geom.Line(A,B), rho);
			Geom.Point Q = Geom.Intersection (new Geom.Line(C,D), rho);
			if (P == null || Q == null)
				return Double.PositiveInfinity; // a special case where one of the SSEs is (almost) parallel to the plane
			return Geom.Distance (P, Q);
		}

		public enum GeometryCheckResult {OK, NOK, IncompleteChain};

		private static bool IsAltLocOK(Atom a){
			return a.AltLoc == Cif.Tables.AtomTable.DEFAULT_ALT_LOC || a.AltLoc == "A";
		}

		public static GeometryCheckResult CheckGeometry (IEnumerable<Residue> residues, char sseType, double RMSDLimit, out double maxRMSD){
			maxRMSD = 0;
			//Lib.WriteDebug ("Geometry check on {0} {1,3} - {2,3}: ", residues.First ().ChainID, residues.First ().ResSeq, residues.Last ().ResSeq);

			List<Atom> CAlphas = residues.SelectMany (r => r.GetCAlphas()).ToList ();
			//Lib.WriteLineDebug ("residues: {0}, CAlphas: {1}", residues.Count(), CAlphas.Count);
			if (CAlphas.Any (a => !IsAltLocOK(a))) {
				if (!alternativeLocationWarningPrinted) {
					Lib.WriteWarning ("Alternative locations found. Ignoring alternative locations except for ' ' and 'A'.");
					alternativeLocationWarningPrinted = true;
				}
				CAlphas = CAlphas.Where (a => IsAltLocOK(a)).ToList ();
			}

			for (int i = 0; i <= CAlphas.Count - 4; i++) {
				List<Atom> quad = CAlphas.GetRange (i, 4);
				//Lib.WriteLineDebug ("quad [{0}]: {1}", quad.Count(), quad.EnumerateWithCommas());
				if (quad.Select ((x, j) => x.ResidueSeqNumber - j).Distinct ().Count () != 1) {
					Lib.WriteLineDebug ("Missing residues somewhere between {0} and {1}", quad [0].ResidueSeqNumber, quad [3].ResidueSeqNumber);
					return GeometryCheckResult.IncompleteChain;
				}
				Matrix mobile = Matrix.FromRowVectors (quad.Select (x => x.Position ()).ToList ());
				Matrix trans;
				Matrix rot;
				double rmsd;
				LibAlgebra.Align (mobile, IdealShapes.GetShape (sseType).Points, out trans, out rot, out rmsd);
				maxRMSD = Math.Max (maxRMSD, rmsd);
				//Console.WriteLine ("<{0,3}> {1,8}", quad [1].ResSeq, rmsd);
				if (rmsd > RMSDLimit) {
					//Lib.WriteLineDebug ("High RMSD ({0}) on {1,3} - {2,3}.", rmsd.ToString ("0.00"), quad [0].ResSeq, quad [3].ResSeq);
					return GeometryCheckResult.NOK;
				}
			}
			//Lib.WriteLineDebug ("OK. (Max RMSD: {0})", maxRMSD.ToString ("0.00"));
			return GeometryCheckResult.OK;
		}

		public static GeometryCheckResult CheckGeometryOf1Unit (IEnumerable<Residue> residues, char sseType, double RMSDLimit, out double rmsd){
			rmsd = 0;
			int unitLength = residues.Count();
			List<Atom> unit = residues
				.Select (r => r.GetCAlphas())
				.TakeWhile(la => la.Count()>0)
				.Select(la => la.First())
				.ToList ();
				
			if (unit.Any (a => !IsAltLocOK(a))) {
				if (!alternativeLocationWarningPrinted) {
					Lib.WriteWarning ("Alternative locations found. Ignoring alternative locations except for ' ' and 'A'.");
					alternativeLocationWarningPrinted = true;
				}
				unit = unit.Where (a => IsAltLocOK(a)).ToList ();
			}
			if (unit.Count != unitLength) {
				return GeometryCheckResult.NOK; //some of the input residues miss CAlpha atom
			}

			if (unit.Select ((x, j) => x.ResidueSeqNumber - j).Distinct ().Count () != 1) {
				Lib.WriteLineDebug ("Missing residues somewhere between {0} and {1}", unit [0].ResidueSeqNumber, unit [3].ResidueSeqNumber);
				return GeometryCheckResult.IncompleteChain;
			}
			Matrix mobile = Matrix.FromRowVectors (unit.Select (x => x.Position ()).ToList ());
			Matrix trans;
			Matrix rot;
			LibAlgebra.Align (mobile, IdealShapes.GetShape (sseType).Points, out trans, out rot, out rmsd);

			return (rmsd <= RMSDLimit) ? GeometryCheckResult.OK : GeometryCheckResult.NOK;
		}

		public static GeometryCheckResult CheckSSEGeometry (Protein p, SSE sse, double RMSDLimit, out double maxRMSD){
			Tuple<int,int>[] ranges = new Tuple<int,int>[]{ new Tuple<int, int> (sse.Start, sse.End) };
			IEnumerable<Residue> residues = p.GetChain (sse.ChainID).GetResidues (ranges).First ();
			return CheckGeometry (residues,sse.Type,RMSDLimit,out maxRMSD);
		}

		public static Vector DirectionVector(IList<Atom> helix){
			if (helix.Count < 3)
				throw new ArgumentException ("Backbone of the helix must contain at least 3 atoms.");
			Vector v = new Vector (0, 0, 0);
			for (int i = 0; i<helix.Count-2; i++) {
				Vector a = helix [i].Position ();
				Vector b = helix [i + 1].Position ();
				Vector c = helix [i + 2].Position ();
				v += ((b - a) % (c - b));
			}
			return -v.Normalize ();
		}


		public static Vector DirectionVector_LinRegVersionAlphaHelix(IList<Atom> helix){
			double rotationPerAtom = 2 * Math.PI / 3.6 / 3;
			if (helix.Count < 3)
				throw new ArgumentException ("Backbone of the helix must contain at least 3 atoms.");
			List<int> indices = helix.ToList().Select ((x, i) => i).ToList();
			List<double> values = indices.Select (i => Math.Cos (rotationPerAtom * i)).ToList ();
			values.AddRange (indices.Select (i => Math.Sin (rotationPerAtom * i)));
			values.AddRange (indices.Select (i => Convert.ToDouble (i)));
			Matrix niceHelixTransposed = Matrix.CreateByRows (3, indices.Count, values);
			niceHelixTransposed.MeanCenterHorizontally ();
			Matrix myHelix = Matrix.FromRowVectors (helix.ToList().Select(x => x.Position()).ToList());
			myHelix.MeanCenterVertically ();
			Matrix rotMatrix = niceHelixTransposed * myHelix;
			rotMatrix.NormalizeRows ();
			if (helix [0].ResidueSeqNumber == 172) {
				Console.WriteLine (rotMatrix.ToRowVectors () [0]);
				Console.WriteLine (rotMatrix.ToRowVectors () [1]);
				Console.WriteLine (rotMatrix.ToRowVectors () [2]);
			}
			return rotMatrix.ToRowVectors () [2];
		}

		/**Calculate a line segment that best describes an SSE given as a list of residues. 
		 * numExtraResidues - indicates the number of residues at both the beginning and the end of the list which are not the part of the SSE but should be used in geometrical calculations.
		 */
		//private static Tuple<Vector,Vector> SSEAsLineSegment_GeomVersion(IEnumerable<Residue> residues, char sseType, int numExtraResidues, out List<double> rmsdList){
		private static Tuple<Vector,Vector> SSEAsLineSegment_GeomVersion(IEnumerable<Residue> residues, SSE sse, out List<double> rmsdList){
			rmsdList = new List<double> ();
				
			if (sse.IsNotFound ()) {
				return new Tuple<Vector,Vector> (Vector.ZERO, Vector.ZERO);
			}

			int startIndex = 0;
			int endIndex = residues.Count () - 1;
			try{
				startIndex = residues.IndicesWhere (r => r.SeqNumber == sse.Start).First ();
				endIndex = residues.IndicesWhere (r => r.SeqNumber == sse.End).First ();
			} catch(Exception) {
				Lib.WriteErrorAndExit ("SSE starting or ending residue is missing: {0}", sse);
			}
			if (residues.Count() < 4 || !residues.All (r => r.HasCAlpha())) {
				Lib.WriteWarning ("At least 4 residues are needed to calculate line segment. Less than 4 residues are given, so atom coordinates will be used instead of line segment approximation.");
				Vector? u_ = residues.ElementAt (startIndex).GetCAlpha()?.Position ();
				Vector? v_ = residues.ElementAt (endIndex).GetCAlpha()?.Position ();
				if (u_ != null && v_ !=null){
					return new Tuple<Vector,Vector> ((Vector) u_, (Vector) v_);
				} else {
					throw new Exception("Important residues are missing C alpha atom.");
				}
			}
			if (residues.Select ((r, i) => r.SeqNumber - i).Distinct ().Count () != 1) {
				Lib.WriteWarning ("Missing residues around " +sse.ToString () + ". Atom coordinates will be used instead of line segment approximation.");
				Vector? u_ = residues.ElementAt (startIndex).GetCAlpha()?.Position ();
				Vector? v_ = residues.ElementAt (endIndex).GetCAlpha()?.Position ();
				if (u_ != null && v_ !=null){
					return new Tuple<Vector,Vector> ((Vector) u_, (Vector) v_);
				} else {
					throw new Exception("Important residues are missing C alpha atom.");
				}
			}

			List<Vector> coords = residues.Select (r => r.GetCAlpha()).Where(a => a != null).Select(a => ((Atom) a).Position()).ToList ();


			Vector sumAxes = Vector.ZERO;
			Vector sumOrigins = Vector.ZERO;
			for (int i = 0; i < coords.Count - 3; i++) {
				Matrix mobile = Matrix.FromRowVectors (coords.GetRange (i, 4));
				Matrix trans;
				Matrix rot;
				double rmsd;
				LibAlgebra.Align (mobile, IdealShapes.GetShape (sse.Type).Points, out trans, out rot, out rmsd);
				Matrix rotBack = rot.Transpose ();
				sumAxes += (IdealShapes.GetShape (sse.Type).Axis * rotBack).ToRowVectors () [0];
				sumOrigins += (IdealShapes.GetShape (sse.Type).Origin * rotBack - trans).ToRowVectors () [0];
				rmsdList.Add (rmsd);
				// if (Double.IsNaN(sumAxes.X)) Lib.WriteLineDebug($"    SumAxes {sumAxes}, SumOrigins {sumOrigins}, Mobile {mobile}, Rot {rot}, Trans {trans}");
				//Console.WriteLine ("<{0,3}> {1,8}", quad [1].ResSeq, rmsd);
			}
			Vector axis = sumAxes.Normalize ();
			Vector origin = sumOrigins / (coords.Count - 3);
			Vector u = origin + ((coords [startIndex] - origin) * axis) * axis;
			Vector v = origin + ((coords [endIndex] - origin) * axis) * axis;
			// Lib.WriteLineDebug($"Axis {axis}, SumAxes {sumAxes}, Origin {origin}, {Double.IsNaN(u.X)}");
			return new Tuple<Vector, Vector> (u, v);
		}

		public static List<SSEInSpace> SSEsAsLineSegments_GeomVersion(Chain chain, List<SSE> sses, out List<double>[] rmsdLists){
			if (sses.Any (x => x.ChainID != chain.Id))
				throw new ArgumentException (System.Reflection.MethodBase.GetCurrentMethod ().Name + ": Some SSEs are not from this chain (" + chain.Id + ").");
			IEnumerable<Tuple<int,int>> ranges = sses.Select (sse => new Tuple<int,int> (sse.Start - NumExtraResidues (sse), sse.End + NumExtraResidues (sse)));
			List<List<Residue>> residueLists = chain.GetResidues (ranges).ToList ();
			for (int i = 0; i < sses.Count; i++) {
				if (residueLists [i].Count != sses [i].Length () + 2 * NumExtraResidues (sses [i])) {
					Lib.WriteWarning (System.Reflection.MethodBase.GetCurrentMethod ().Name + ": Some residues might be missing between " + (sses [i].Start-NumExtraResidues (sses[i])) + " and " + (sses [i].End+NumExtraResidues (sses[i])) + " in chain " + sses [i].ChainID + ": " + residueLists [i].Aggregate ("", (s, r) => s + ", " + r.ToString (), s => s.Length>=2 ? s.Substring (2):s));
				}
			} 
			//for (int i = 0; i < residueLists.Count; i++) { Console.WriteLine ("{0} in {1} {2}-{3}: {4}", sses [i].Label, sses [i].ChainID, sses [i].Start, sses [i].End, residueLists [i].Aggregate ("", (x, y) => x + ", " + y.ToString ())); }
			List<double>[] rmsdLists_ = new List<double>[sses.Count];
			List<SSEInSpace> ssesInSpace = 
				sses.Select ((sse, i) => new SSEInSpace (sse, LibAnnotation.SSEAsLineSegment_GeomVersion (residueLists [i], sse, out rmsdLists_ [i]))).ToList ();
			rmsdLists = rmsdLists_;
			return ssesInSpace;
		}

		private static int NumExtraResidues(SSE sse){
			return (sse.Length () >= 4) ? 0 : (sse.Length () >= 2) ? 1 : 2;
		}

		public static double MetricNo3Neg(SSEInSpace sse1,SSEInSpace sse2){
			//this metric is to be maximized (minimized in absolute value)
			return -MetricNo3Pos (sse1, sse2);
		}

		public static double MetricNo3Pos(SSEInSpace sse1,SSEInSpace sse2){
			//this metric is to be minimized
			return (sse1.StartVector - sse2.StartVector).Size 
				+ (sse1.EndVector - sse2.EndVector).Size;
		}

		public static double MetricNo4(SSEInSpace sse1,SSEInSpace sse2){
			//this metric is to be maximized (minimized in absolute value)
			return - (sse1.StartVector-sse2.StartVector).SqSize
				- (sse1.StartVector-sse2.StartVector).SqSize;
		}

		public static double MetricNo5Pos(SSEInSpace sse1,SSEInSpace sse2){
			//this metric is to be minimized
			double m3 = (sse1.StartVector - sse2.StartVector).Size
			            + (sse1.EndVector - sse2.EndVector).Size;
			Vector u = sse1.EndVector - sse1.StartVector;
			Vector v = sse2.EndVector - sse2.StartVector;
			double angleScaling = 2.0 - u * v / u.Size / v.Size;
			return m3 * angleScaling;
		}

		private static double PartialMetricNo6(Vector[] ra, Vector[] rb, int d){
			double sum = 0;
			//extras in ra before
			for (int i = 0; i < d; i++)
				sum += (ra [i] - rb [0]).Size;
			//extras in rb before
			for (int i = d; i < 0; i++)
				sum += (ra [0] - rb [i - d]).Size;
			//overlapping part
			for (int i = Math.Max (0, d); i < Math.Min (ra.Length, d + rb.Length); i++)
				sum += (ra [i] - rb [i - d]).Size;
			//extras in ra after
			for (int i = d+rb.Length; i < ra.Length; i++)
				sum += (ra [i] - rb [rb.Length-1]).Size;
			//extras in rb after
			for (int i = ra.Length; i < d+rb.Length; i++)
				sum += (ra [ra.Length-1] - rb [i - d]).Size;
			//Lib.WriteDebug (Math.Max (ra.Length, rb.Length + d) - Math.Min (0, d)+", ");
			return sum / (Math.Max (ra.Length, rb.Length + d) - Math.Min (0, d));
		}

		public static double MetricNo6Pos(SSE sse1, SSE sse2, Chain chain1, Chain chain2){
			Vector[] ra = chain1.GetResidues (new Tuple<int,int>[] { new Tuple<int, int> (sse1.Start, sse1.End) })[0].Select (r=>r.GetCAlphas().First().Position ()).ToArray ();
			Vector[] rb = chain2.GetResidues (new Tuple<int,int>[] { new Tuple<int, int> (sse2.Start, sse2.End) })[0].Select (r=>r.GetCAlphas().First().Position ()).ToArray ();
			int minDisplacement = -rb.Length + 1;
			int maxDisplacement = ra.Length - 1;
			IEnumerable<int> displacements = Enumerable.Range (minDisplacement, maxDisplacement-minDisplacement+1);
			//Lib.WriteLineDebug ("\n{0} [{1}] {2} [{3}] ",sse1, ra.Length,sse2, rb.Length);
			//Lib.WriteLineDebug ("Displacements: {0}", Lib.EnumerateWithCommas (displacements));
			return displacements.Select (d => PartialMetricNo6 (ra, rb, d)).Min ();
		}

		public static double MetricNo7Pos(SSE sse1, SSE sse2, List<Tuple<int?,int?,double>> alignment, Dictionary<int,int> tResiToAli, Dictionary<int,int> qResiToAli){
			try{
			int s1 = tResiToAli [sse1.Start];
			int s2 = qResiToAli [sse2.Start];
			int e1 = tResiToAli [sse1.End];
			int e2 = qResiToAli [sse2.End];
			int deltaStart = Math.Abs (s1-s2);
			if (s1 < s2) {
				int j = s2;
				while (j < alignment.Count && alignment [j].Item1 == null)
					j++;
				if (j < alignment.Count)
					deltaStart = Math.Min (deltaStart, (int)alignment [j].Item1 - (int)alignment [s1].Item1);
			} else {
				int j = s1;
				while (j < alignment.Count && alignment [j].Item2 == null)
					j++;
				if (j < alignment.Count)
					deltaStart = Math.Min (deltaStart, (int)alignment [j].Item2 - (int)alignment [s2].Item2);
			}
			int deltaEnd = Math.Abs (e1-e2);
			if (e1 > e2) {
				int j = e2;
				while (j >= 0 && alignment [j].Item1 == null)
					j--;
				if (j >= 0)
					deltaEnd = Math.Min (deltaEnd, (int)alignment [e1].Item1 - (int)alignment [j].Item1);
			} else {
				int j = e1;
				while (j >= 0 && alignment [j].Item2 == null)
					j--;
				if (j >= 0)
					deltaEnd = Math.Min (deltaEnd, (int)alignment [e2].Item2 - (int)alignment [j].Item2);
			}

			return deltaStart + deltaEnd;
			}catch(Exception e){
				Lib.WriteLineDebug ("m7: {0} {1} {2} {3}", sse1.Start, sse2.Start, sse1.End, sse2.End);
				throw e;
			}
		}

		public static double LengthDiffPenalty(SSE sse1, SSE sse2){
			const double ALPHA = 10;
			const double BETA_SQ = 9;
			double len1 = sse1.End - sse1.Start + 1;
			double len2 = sse2.End - sse2.Start + 1;
			return ALPHA * Math.Abs (len1 - len2) / Math.Sqrt (len1*len2 + BETA_SQ);
		}

		/*public static int FindBestCounterpart(Tuple<Vector,Vector> myHelix, 
			List<Tuple<Vector,Vector>> candidateHelices, 
			Func<Tuple<Vector,Vector>,Tuple<Vector,Vector>,double> metric)
		{
			if (candidateHelices.Count < 1)
				throw new InvalidDataException ("At least 1 candidate helix needed.");
			int best=0;
			double max=Double.MinValue;
			for (int i=0; i<candidateHelices.Count; i++) {	
				double value = metric (myHelix, candidateHelices [i]);
				if (value > max) {
					best = i;
					max = value;
				}
			}
			return best;
		}*/

		private class SSET<T_>{
			public SSE SSE{ get; set; }
			public T_ T { get; set; }
			public SSET(SSE sse, T_ t){
				SSE=sse;
				T=t;
			}
		}

		private class SSETint<T_>{
			public SSE SSE{ get; set; }
			public T_ T { get; set; }
			public int Int { get; set; }
			public SSETint(SSE sse, T_ t, int i){
				SSE=sse;
				T=t;
				Int=i;
			}
		}

		private class SSEInSpaceInt{
			public SSEInSpace SSE{ get; set; }
			public int Int { get; set; }
			public SSEInSpaceInt(SSEInSpace sse, int i){
				SSE=sse;
				Int=i;
			}
		}

		/** Returns index of best counterpart in candidateHelices and the value of metric. Returns null if not found. */
		private static Tuple<int,double> FindBestCounterpartNew(SSEInSpaceInt myHelixT,List<SSEInSpaceInt> candidateHelices, int iFrom, int iTo, Func<char,char,bool> typeMatching, Func<SSEInSpace,SSEInSpace,double> metricToMax)
		{
			if (iFrom > iTo || iFrom > candidateHelices.Count - 1 || iTo < 0)
				return null; // not found
			int best=-1; // not found
			double max=Double.MinValue;
			for (int i=iFrom; i<=iTo; i++) {	
				if (candidateHelices [i].Int == -1 && typeMatching (myHelixT.SSE.Type, candidateHelices [i].SSE.Type)) { // not annotated yet and having matching type
					double value = metricToMax (myHelixT.SSE, candidateHelices [i].SSE);
					if (value > max) { // so far the best
						best = i;
						max = value;
					}
				}
			}
			return best!=-1 ? new Tuple<int, double>(best,max) : null;
		}

		/** For each SSE in templates finds the best-fitting SSE in candidates. 
			templates and candidates bear additional information (of type T) used for calculation of the metric.
			chainMapping indicates SSEs from which template chain should be looked for in which protein chain.
			SSEs are choose in such way that the value of metric is minimized.
			*/
		public static List<Tuple<SSEInSpace,double>> FindCounterpartsNew(IEnumerable<SSEInSpace> templates, IEnumerable<SSEInSpace> candidates, Dictionary<string,string> chainMapping, Func<char,char,bool> typeMatching, Func<SSEInSpace,SSEInSpace,double> metricToMax){

			Dictionary<string,List<SSEInSpaceInt>> templatesByResseq = new Dictionary<string,List<SSEInSpaceInt>> (); // list of SSEs for each template chain sorted by starting resseq
			Dictionary<string,List<SSEInSpaceInt>> templatesByPriority = new Dictionary<string,List<SSEInSpaceInt>> (); // list of SSEs for each template chain sorted by their order of processing
			Dictionary<string,List<SSEInSpaceInt>> candidatesByResseq = new Dictionary<string,List<SSEInSpaceInt>> (); // list of SSEs for each template chain sorted by starting resseq
			List<Tuple<SSEInSpace,double>> result = new List<Tuple<SSEInSpace,double>> ();

			foreach (string chain in chainMapping.Keys) {
				templatesByResseq [chain] = new List<SSEInSpaceInt> ();
				templatesByPriority [chain] = new List<SSEInSpaceInt> ();
				candidatesByResseq [chainMapping [chain]] = new List<SSEInSpaceInt> ();
			}

			foreach (SSEInSpace template in templates)
				if (templatesByResseq.ContainsKey (template.ChainID)) {
					templatesByResseq [template.ChainID].Add (new SSEInSpaceInt (template, -1)); //-1=no counterpart selected
					templatesByPriority [template.ChainID].Add (new SSEInSpaceInt (template, -1)); //-1=undefined
				}

			foreach (SSEInSpace candidate in candidates)
				if (candidatesByResseq.ContainsKey (candidate.ChainID))
					candidatesByResseq [candidate.ChainID].Add (new SSEInSpaceInt (candidate, -1)); // -1=not annotated yet

			foreach (string chain in chainMapping.Keys) {
				templatesByResseq [chain].Sort ((x, y) => x.SSE.Start.CompareTo (y.SSE.Start));
				candidatesByResseq [chainMapping [chain]].Sort ((x, y) => x.SSE.Start.CompareTo (y.SSE.Start));
				/*templatesByResseq [chain].ForEach (x => Console.WriteLine (x.SSE.ToString ()));
				candidatesByResseq [chainMapping[chain]].ForEach (x => Console.WriteLine (x.SSE.ToString ()));*/
			}

			foreach (SSEInSpace t in templates.Where(x=>chainMapping.ContainsKey(x.ChainID))) {
				var temps = templatesByResseq [t.ChainID];
				var cands = candidatesByResseq [chainMapping [t.ChainID]];
				int i = temps.FindIndex (x => x.SSE.Equals(t));
				int before = i;
				int after = i;
				while (before >= 0 && temps [before].Int == -1)
					before--;
				while (after < temps.Count && temps [after].Int == -1)
					after++;
				/*Console.WriteLine ("Searching " + temps [i].SSE.Label);
				cands.ForEach (x => Console.Write (x.SSE.Start.ToString ()+ "\t"));
				Console.WriteLine ();
				cands.ForEach (x => Console.Write (before>=0&&cands.IndexOf(x)==temps[before].Int+1?">\t":after<temps.Count&&cands.IndexOf(x)==temps[after].Int-1?"<\t":x.Int==-1?"\t":x.Int.ToString () + "\t"));
				Console.WriteLine ();
				cands.ForEach (x => Console.Write (metric(x.SSE,t).ToString()+ "\t"));
				Console.WriteLine ();*/
				Tuple<int,double> match = FindBestCounterpartNew (temps [i], cands, before>=0 ? temps[before].Int + 1 : 0, after<temps.Count ? temps[after].Int - 1 : cands.Count-1, typeMatching, metricToMax);
				if (match!=null) {
					temps [i].Int = match.Item1;
					cands [match.Item1].Int = i; // already annotated
					result.Add(new Tuple<SSEInSpace,double>(cands[match.Item1].SSE.RelabeledCopy(t.Label),match.Item2));
				} else {
					result.Add (new Tuple<SSEInSpace,double>(SSEInSpace.NewNotFound(t.Label), Double.NegativeInfinity)); // counterpart not found
				}
				/*Console.WriteLine ("Found at " + (match==null? "nowhere": match.Item1.ToString())+".");
				Console.WriteLine ();*/
			}

			return result;
		}

		/*
		private static Tuple<int,double> FindBestCounterpartNew<T>(T myHelixT,List<SSETint<T>> candidateHelices, int iFrom, int iTo, Func<T,T,double> metric)
		{
			if (iFrom > iTo || iFrom > candidateHelices.Count - 1 || iTo < 0)
				return null; // not found
			int best=-1; // not found
			double max=Double.MinValue;
			for (int i=iFrom; i<=iTo; i++) {	
				double value = metric (myHelixT, candidateHelices [i].T);
				if (candidateHelices[i].Int==-1 && value > max) { // so far the best and not annotated yet
					best = i;
					max = value;
				}
			}
			return best!=-1 ? new Tuple<int, double>(best,max) : null;
		}
		public static List<Tuple<SSE,double>> FindCounterpartsNew<T>(IEnumerable<Tuple<SSE,T>> templates, IEnumerable<Tuple<SSE,T>>candidates, Dictionary<char,char> chainMapping, Func<T,T,double> metric){
			//TODO Test

			Dictionary<char,List<SSETint<T>>> templatesByResseq = new Dictionary<char,List<SSETint<T>>> (); // list of SSEs for each template chain sorted by starting resseq
			Dictionary<char,List<SSETint<T>>> templatesByPriority = new Dictionary<char,List<SSETint<T>>> (); // list of SSEs for each template chain sorted by their order of processing
			Dictionary<char,List<SSETint<T>>> candidatesByResseq = new Dictionary<char,List<SSETint<T>>> (); // list of SSEs for each template chain sorted by starting resseq
			List<Tuple<SSE,double>> result = new List<Tuple<SSE,double>> ();

			foreach (char chain in chainMapping.Keys) {
				templatesByResseq [chain] = new List<SSETint<T>> ();
				templatesByPriority [chain] = new List<SSETint<T>> ();
				candidatesByResseq [chainMapping [chain]] = new List<SSETint<T>> ();
			}

			foreach (Tuple<SSE,T> template in templates)
				if (templatesByResseq.ContainsKey (template.Item1.ChainID)) {
					templatesByResseq [template.Item1.ChainID].Add (new SSETint<T> (template.Item1, template.Item2, -1)); //-1=no counterpart selected
					templatesByPriority [template.Item1.ChainID].Add (new SSETint<T> (template.Item1, template.Item2, -1)); //-1=undefined
				}

			foreach (Tuple<SSE,T> candidate in candidates)
				if (candidatesByResseq.ContainsKey (candidate.Item1.ChainID))
					candidatesByResseq [candidate.Item1.ChainID].Add (new SSETint<T> (candidate.Item1, candidate.Item2, -1)); // -1=not annotated yet

			foreach (char chain in chainMapping.Keys) {
				templatesByResseq [chain].Sort ((x, y) => x.SSE.Start.CompareTo (y.SSE.Start));
				candidatesByResseq [chainMapping [chain]].Sort ((x, y) => x.SSE.Start.CompareTo (y.SSE.Start));
				//templatesByResseq [chain].ForEach (x => Console.WriteLine (x.SSE.ToString ()));
				//candidatesByResseq [chainMapping[chain]].ForEach (x => Console.WriteLine (x.SSE.ToString ()));
			}

			foreach (SSE t in templates.Select(x=>x.Item1).Where(x=>chainMapping.ContainsKey(x.ChainID))) {
				var temps = templatesByResseq [t.ChainID];
				var cands = candidatesByResseq [chainMapping [t.ChainID]];
				int i = temps.FindIndex (x => x.SSE.Equals(t));
				int before = i;
				int after = i;
				while (before >= 0 && temps [before].Int == -1)
					before--;
				while (after < temps.Count && temps [after].Int == -1)
					after++;
				//temps.ForEach(x=>Console.Write ("{0}, ", x.Int));
				//Console.WriteLine ("\ni = {0}, before = {1}, after = {2}, temps.Count = {3} \n", i,before,after,temps.Count);
				Tuple<int,double> match = FindBestCounterpartNew (temps [i].T, cands, before>=0 ? temps[before].Int + 1 : 0, after<temps.Count ? temps[after].Int - 1 : cands.Count-1, metric);
				if (match!=null) {
					temps [i].Int = match.Item1;
					cands [match.Item1].Int = i; // already annotated
					result.Add (new Tuple<SSE,double>(new SSE (t.Label, cands [match.Item1].SSE.ChainID, cands [match.Item1].SSE.Start, cands [match.Item1].SSE.End, cands [match.Item1].SSE.Type), match.Item2));
				} else {
					result.Add (new Tuple<SSE,double>(new SSE (t.Label, chainMapping[t.ChainID], NotFoundStart, NotFoundEnd, NotFoundType), Double.NegativeInfinity)); // counterpart not found
				}
			}

			return result;
		}*/

		public static List<Tuple<SSEInSpace,double>> FindCounterparts_DynProg(IList<SSEInSpace> templates, IList<SSEInSpace> candidates,  Func<char,char,bool> typeMatching, 
			Func<SSEInSpace,SSEInSpace,double> metricToMin,Func<SSEInSpace,double> skipTemplatePenalty,Func<SSEInSpace,double> skipCandidatePenalty, out List<double> suspiciousnessList){

			if (templates.Select (x => x.ChainID).Distinct ().Count () > 1)
				Lib.WriteWarning ("{0}: passed template SSEs from more than one chain.", System.Reflection.MethodBase.GetCurrentMethod ().Name);
			if (candidates.Select (x => x.ChainID).Distinct ().Count () > 1)
				Lib.WriteWarning ("{0}: passed candidate SSEs from more than one chain.", System.Reflection.MethodBase.GetCurrentMethod ().Name);

			int[] templateRanks;
			SSEInSpace[] templateArray = Lib.OrderAndGetRanks (templates, sse => sse.Start, out templateRanks).ToArray ();
			int[] candidateRanks;
			SSEInSpace[] candidateArray = Lib.OrderAndGetRanks (candidates, sse => sse.Start, out candidateRanks).ToArray ();
			SSEInSpace[] selectedCandidateArray = new SSEInSpace[templateArray.Length];
			int[] selectedJ = new int[templateArray.Length];
			double[] metricArray = new double[templateArray.Length];

			int m = templateArray.Length + 1;
			int n = candidateArray.Length + 1;


			double[,] metricMatrix = new double[m-1, n-1];
			double[,] dynProgMatrix = new double[m, n];
			int[,] reconstrMatrix = new int[m, n];
			//coding reconstrMatrix: 
			//  -1: optimum was obtained from the left neighbor (skip candidate SSE)
			//   1: optimum was obtained from upper neighbor (skip template SSE)
			//   0: optimum was obtained from diagonal neighbor (pair a template and candidate SSE)

			// Calculation of the matrix of metric.
			for (int i = 0; i < m-1; i++) {
				for (int j = 0; j < n-1; j++) {
					metricMatrix [i, j] = metricToMin (templateArray [i], candidateArray [j]);
				}
			}

			// Calculation of dynamic programming matrix.
			dynProgMatrix [0, 0] = 0.0;

			for (int i = 1; i < m; i++) {
				reconstrMatrix [i, 0] = 1;
				dynProgMatrix [i, 0] = dynProgMatrix [i - 1, 0] + skipTemplatePenalty (templateArray [i - 1]);
			}

			for (int j = 1; j < n; j++) {
				reconstrMatrix [0, j] = -1;
				dynProgMatrix [0, j] = dynProgMatrix [0, j - 1] + skipCandidatePenalty (candidateArray [j - 1]);
			}

			for (int i = 1; i < m; i++) {
				for (int j = 1; j < n; j++) {
					SSEInSpace temp = templateArray [i - 1];
					SSEInSpace cand = candidateArray [j - 1];

					double valueUp = dynProgMatrix [i - 1, j] + skipTemplatePenalty (temp);
					double valueLeft = dynProgMatrix [i, j - 1] + skipCandidatePenalty (cand);

					if (valueUp <= valueLeft) {
						//skip template
						reconstrMatrix [i, j] = 1;
						dynProgMatrix [i, j] = valueUp;
					} else {
						//skip candidate
						reconstrMatrix [i, j] = -1;
						dynProgMatrix [i, j] = valueLeft;
					}

					if (typeMatching (temp.Type, cand.Type)) {
						double valueDiag = dynProgMatrix [i - 1, j - 1] + metricToMin (temp, cand);
						if (valueDiag <= dynProgMatrix [i, j]) {
							//pair template with candidate
							reconstrMatrix [i, j] = 0;
							dynProgMatrix [i, j] = valueDiag;
						}
					}
				}
			}


			// print raw metric matrix 
			if (Lib.DoWriteDebug) {
				TextWriter w = new StreamWriter (Path.Combine (MainClass.Directory, "metric_matrix.tsv"));
				w.Write ("\t<length>\t");
				for (int j = 0; j < candidateArray.Length; j++) {
					w.Write (candidateArray [j].Label + "\t");
				}
				w.WriteLine ();
				for (int i = 0; i < templateArray.Length; i++) {
					w.Write (templateArray [i].Label + "\t");
					w.Write ((templateArray [i].EndVector - templateArray [i].StartVector).Size.ToString ("0.0") + "\t");
					for (int j = 0; j < candidateArray.Length; j++) {
						w.Write (metricToMin (templateArray [i], candidateArray [j]).ToString ("0.0") + "\t");
					}
					w.WriteLine ();
				}
				w.WriteLine ();
				w.WriteLine ();
				w.Close ();
			}

			// print dynamic programming matrix
			if (Lib.DoWriteDebug) {
				TextWriter w = new StreamWriter (Path.Combine (MainClass.Directory, "dyn_prog_matrix.tsv"));
				w.Write ("\t");
				for (int j = 0; j < n; j++) {
					w.Write (j == 0 ? "-\t" : candidateArray [j - 1].Label + "\t");
				}
				w.WriteLine ();
				for (int i = 0; i < m; i++) {
					w.Write ("\t");
					for (int j = 0; j < n; j++) {
						w.Write (reconstrMatrix [i, j] == -1 ? ">" : reconstrMatrix [i, j] == 1 ? "v" : "\\");
						w.Write ("\t");
					}
					w.WriteLine ();
					w.Write (i == 0 ? "-\t" : templateArray [i - 1].Label + "\t");
					for (int j = 0; j < n; j++) {
						w.Write (dynProgMatrix [i, j].ToString ("0.0") + "\t");
					}
					w.WriteLine ();
				}
				w.Close ();
			}

			// Reconstruction of the best solution.
			int row = m - 1;
			int col = n - 1;
			while (row != 0 || col != 0) {
				switch (reconstrMatrix [row, col]) {
				case -1:
					// skipped candidate SSE
					col--;
					break;
				case 1:
					// skipped template SSE
					selectedCandidateArray [row - 1] = SSEInSpace.NewNotFound (templateArray [row - 1].Label);
					selectedJ [row - 1] = -1;
					metricArray [row - 1] = Double.PositiveInfinity;
					row--;
					break;
				case 0:
					// paired template and candidate SSE
					selectedCandidateArray [row - 1] = candidateArray [col - 1].RelabeledCopy (templateArray [row - 1].Label);
					selectedJ [row - 1] = col - 1;
					metricArray [row - 1] = dynProgMatrix [row, col] - dynProgMatrix [row - 1, col - 1];
					row--;
					col--;
					break;
				}
			}

			// Finding possibly ambiguous annotations.
			double[] suspiciousnessArray=new double[m-1];
			for(int i=0;i<m-1;i++) {
				if (selectedJ [i] >= 0) {
					double otherMinInRow = Enumerable.Range (0, n - 1)
						.Where (j => typeMatching (templateArray [i].Type, candidateArray [j].Type))
						.Where (j => j != selectedJ [i])
						.Select (j => metricMatrix [i, j])
						.Min ();
					double otherMinInColumn = Enumerable.Range (0, m - 1)
						.Where (j => typeMatching (templateArray [j].Type, candidateArray [selectedJ [i]].Type))
						.Where (j => j != i)
						.Select (j => metricMatrix [j, selectedJ [i]])
						.Min ();
					double otherMin = Math.Min (otherMinInRow, otherMinInColumn);
					suspiciousnessArray [i] = metricArray [i] / (otherMin + metricArray [i]);
					if (suspiciousnessArray [i] >= SUSPICIOUSNESS_THRESHOLD) {
						Lib.WriteWarning ("Possibly ambiguous annotation around \"{0}\" (annotated / (alternative+annotated) = {1}).", templateArray [i].Label, suspiciousnessArray [i].ToString ("0.00"));
					}
				} else {
					suspiciousnessArray [i] = Double.NaN;
				}
			}

			//checking sheet ID correspondence
			if (templateArray.Any (sse => sse.IsSheet && sse.SheetId == null) || selectedCandidateArray.Any (sse => sse.IsSheet && sse.SheetId == null)) {
				//skip checking (some strand miss sheet ID)
			} else {
				Dictionary<int?,List<int?>> dictSheetIdByTemplate = new Dictionary<int?,List<int?>> ();
				Dictionary<int?,List<int?>> dictSheetIdByCandidate = new Dictionary<int?,List<int?>> ();
				foreach (int rank in templateRanks) {
					SSE template = templateArray [rank];
					SSE candidate = selectedCandidateArray [rank];
					if (template.IsSheet && candidate.IsSheet) {
						dictSheetIdByTemplate.MultidictionaryAdd (template.SheetId, candidate.SheetId);
						dictSheetIdByCandidate.MultidictionaryAdd (candidate.SheetId, template.SheetId);
					}
				}
				if (dictSheetIdByTemplate.Values.Any (l => l.Distinct ().Count () > 1) || dictSheetIdByCandidate.Values.Any (l => l.Distinct ().Count () > 1)) {
					// sheet ID mismatch found
					Lib.WriteWarning ("Some beta-strands caused a sheet ID mismatch!");
					foreach (var kv in dictSheetIdByTemplate.Where (kv => kv.Value.Distinct ().Count () > 1)) {
						Lib.WriteWarning ("Template sheet {0} matched to query sheets {1}. ", kv.Key, kv.Value.Distinct ().EnumerateWithCommas ());
					}
					foreach (var kv in dictSheetIdByCandidate.Where (kv => kv.Value.Distinct ().Count () > 1)) {
						Lib.WriteWarning ("Query sheet {0} matched to template sheets {1}. ", kv.Key, kv.Value.Distinct ().EnumerateWithCommas ());
					}
				} else {
					// OK -> renaming sheet IDs
					foreach (int rank in templateRanks) {
						SSE template = templateArray [rank];
						SSE candidate = selectedCandidateArray [rank];
						if (template.IsSheet && candidate.IsSheet) {
							candidate.SheetId = dictSheetIdByTemplate [template.SheetId] [0];
						}
					}
				}
			}

			suspiciousnessList = templateRanks.Select (rank => suspiciousnessArray [rank]).ToList ();
			return templateRanks.Select (rank => new Tuple<SSEInSpace,double> (selectedCandidateArray [rank], metricArray [rank])).ToList ();
		}

		/*public static List<Tuple<SSEInSpace,double>> FindCounterparts_BetaGraphMatching(IList<SSEInSpace> templates, IList<SSEInSpace> candidates,
			IEnumerable<Tuple<int,int>> templateNeighbours, IEnumerable<Tuple<int,int>> candidateNeighbours, Func<char,char,bool> typeMatching, 
			Func<SSEInSpace,SSEInSpace,double> metricToMin,Func<SSEInSpace,double> skipTemplatePenalty, Func<SSEInSpace,double> skipCandidatePenalty, out List<double> suspiciousnessList){

			if (templates.Select (x => x.ChainID).Distinct ().Count () > 1)
				Lib.WriteWarning ("{0}: passed template SSEs from more than one chain.", System.Reflection.MethodBase.GetCurrentMethod ().Name);
			if (candidates.Select (x => x.ChainID).Distinct ().Count () > 1)
				Lib.WriteWarning ("{0}: passed candidate SSEs from more than one chain.", System.Reflection.MethodBase.GetCurrentMethod ().Name);
			
			int[] ranksT;
			SSEInSpace[] verticesT = Lib.OrderAndGetRanks (templates, sse => sse.Start, out ranksT).ToArray ();
			int[] ranksQ;
			SSEInSpace[] verticesQ = Lib.OrderAndGetRanks (candidates, sse => sse.Start, out ranksQ).ToArray ();

			int nT = verticesT.Length;
			int nQ = verticesQ.Length;

			bool[,] edgesT=new bool[nT,nT];
			foreach (var neighbours in templateNeighbours) {
				edgesT [ranksT [neighbours.Item1], ranksT [neighbours.Item2]] = true;
				edgesT [ranksT [neighbours.Item2], ranksT [neighbours.Item1]] = true;
			}

			bool[,] edgesQ=new bool[nQ,nQ];
			foreach (var neighbours in candidateNeighbours) {
				edgesQ [ranksQ [neighbours.Item1], ranksQ [neighbours.Item2]] = true;
				edgesQ [ranksQ [neighbours.Item2], ranksQ [neighbours.Item1]] = true;
			}

			double[,] scores = new double[nT,nQ];
			for (int i = 0; i < nT; i++) {
				for (int j = 0; j < nQ; j++) {
					SSEInSpace x = verticesT [i];
					SSEInSpace y = verticesQ [j];
					scores [i, j] = typeMatching (x.Type, y.Type) ? 
						skipTemplatePenalty (x) + skipCandidatePenalty (y) - metricToMin (x, y) 
						: 0;
				}
			}

			DateTime stamp = DateTime.Now;
			Lib.WriteLineDebug ("FindCounterparts_BetaGraphMatching: initialized - {0} vs. {1} vertices ({2})",nT,nQ,stamp);

			List<Tuple<int,int>> bestMatching = MaxWeightOrderedMatching (nT, nQ, edgesT, edgesQ, scores);

			Lib.WriteLineDebug ("FindCounterparts_BetaGraphMatching: found matching ({0})",DateTime.Now);
			Lib.WriteLineDebug ("FindCounterparts_BetaGraphMatching: time: {0}",DateTime.Now-stamp);

			SSEInSpace[] selectedVerticesQ = new SSEInSpace[nT];
			foreach (var match in bestMatching) {
				selectedVerticesQ [match.Item1] = verticesQ [match.Item2].RelabeledCopy (verticesT [match.Item1].Label);
			}

			suspiciousnessList = templates.Select (x => Double.NaN).ToList ();
			return ranksT
				.Where (rank=>selectedVerticesQ [rank] != null)
				.Select (rank =>new Tuple<SSEInSpace,double> (selectedVerticesQ [rank], metricToMin (verticesT [rank], selectedVerticesQ [rank])))
				.ToList ();
		}*/


		public static Tuple<List<int>,List<int>> AlignmentToPositions (List<Tuple<int?,int?,double>> alignment){
			List<int> tPositions = new List<int> ();
			List<int> qPositions = new List<int> ();
			for (int i = 0; i < alignment.Count; i++) {
				if (alignment [i].Item1 != null)
					tPositions.Add (i);
				if (alignment [i].Item2 != null)
					qPositions.Add (i);
			}
			return new Tuple<List<int>, List<int>> (tPositions, qPositions);
		}

		public static Dictionary<int,int> DictResiToAli(List<Residue> residues, List<int> positions){
			if (residues.Count != positions.Count)
				throw new ArgumentException ();
			Dictionary<int,int> dictResiToAli = new Dictionary<int,int> ();
			for (int i = 0; i < residues.Count; i++) {
				dictResiToAli.Add (residues[i].SeqNumber,positions[i]);
			}
			return dictResiToAli;
		}

		public static List<Tuple<int?,int?,double>> AlignPoints_DynProg(int m, int n, Func<int,int,double> scoreFunction){
			double[,] scoreMatrix = new double[m, n];
			double[,] dynprogMatrix = new double[m+1, n+1];
			int[,] reconstrMatrix = new int[m+1, n+1];
			List<Tuple<int?,int?,double>> alignment = new List<Tuple<int?,int?,double>>();
			//coding reconstrMatrix: 
			//  -1: optimum was obtained from the left neighbor (skip candidate SSE)
			//   1: optimum was obtained from upper neighbor (skip template SSE)
			//   0: optimum was obtained from diagonal neighbor (pair a template and candidate SSE)

			// Calculation of score matrix.
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) {
					scoreMatrix [i, j] = scoreFunction (i, j);
				}
			}

			// Calculation of dynamic programming matrix.

			for (int i = 0; i <= m; i++) {
				reconstrMatrix [i, 0] = 1;
				dynprogMatrix [i, 0] = 0.0;
			}

			for (int j = 1; j <= n; j++) {
				reconstrMatrix [0, j] = -1;
				dynprogMatrix [0, j] = 0.0;
			}

			for (int i = 1; i <= m; i++) {
				for (int j = 1; j <= n; j++) {
					double valueUp = dynprogMatrix [i - 1, j];
					double valueLeft = dynprogMatrix [i, j - 1];
					double valueDiag = dynprogMatrix [i - 1, j - 1] + scoreMatrix[i-1, j-1];

					if (valueDiag >= valueLeft) {
						if (valueDiag >= valueUp) {
							//diag (match)
							reconstrMatrix [i, j] = 0;
							dynprogMatrix [i, j] = valueDiag;
						} else {
							//up (skip template)
							reconstrMatrix [i, j] = 1;
							dynprogMatrix [i, j] = valueUp;
						}
					} else {
						if (valueLeft >= valueUp) {
							//left (skip candidate)
							reconstrMatrix [i, j] = -1;
							dynprogMatrix [i, j] = valueLeft;
						}else{
							//up (skip template)
							reconstrMatrix [i, j] = 1;
							dynprogMatrix [i, j] = valueUp;
						}
					}
				}
			}

			// print score matrix and dynamic programming matrix
			if (Lib.DoWriteDebug) {
				StreamWriter w = new StreamWriter (Path.Combine (MainClass.Directory, "score_matrix-residue_alignment.tsv"));
				for (int i = 0; i < m; i++) {
					for (int j = 0; j < n; j++) {
						w.Write (scoreMatrix[i, j] + "\t");
					}
					w.WriteLine ();
				}
				w.Close ();
			}
			if (Lib.DoWriteDebug) {
				StreamWriter w = new StreamWriter (Path.Combine (MainClass.Directory, "dynprog_matrix-residue_alignment.tsv"));
				for (int i = 0; i <= m; i++) {
					for (int j = 0; j <= n; j++) {
						w.Write (dynprogMatrix[i,j] + "\t");
					}
					w.WriteLine ();
				}
				w.Close ();
			}

			// Reconstruction of the best solution.
			int row = m;
			int col = n;
			int alignmentPosition = 0;
			while (row != 0 || col != 0) {
				switch (reconstrMatrix [row, col]) {
				case -1:
					// skipped candidate SSE
					alignment.Add (new Tuple<int?, int?, double> (null, col-1, 0));
					alignmentPosition--;
					col--;
					break;
				case 1:
					// skipped template SSE
					alignment.Add (new Tuple<int?, int?, double> (row-1, null,0));
					alignmentPosition--;
					row--;
					break;
				case 0:
					// paired template and candidate SSE
					alignment.Add (new Tuple<int?,int?,double> (row-1, col-1, scoreMatrix[row-1,col-1]));
					alignmentPosition--;
					row--;
					col--;
					break;
				}
			}

			alignment.Reverse ();

			return alignment;
		}

		public static List<Tuple<int?,int?,double>> AlignResidues_DynProg(IList<Residue> templates, IList<Residue> candidates){
			if (templates.Any (r => !r.HasCAlpha()))
				throw new ArgumentException ("Some template residues do not have C-alpha atom.");
			if (candidates.Any (r => !r.HasCAlpha()))
				throw new ArgumentException ("Some candidate residues do not have C-alpha atom.");

			Lib.WriteLineDebug ("Align residues: Candidate residues: {0}", candidates.Select (r => r.ToString ()).EnumerateWithCommas ());
			
			Vector[] tVectors = templates.Select (r => r.GetCAlphas().First().Position ()).ToArray ();
			Vector[] qVectors = candidates.Select (r => r.GetCAlphas().First().Position ()).ToArray ();

			if (Lib.DoWriteDebug) {
				//TODO implement this somehow!
				ModelBuilder builder = new ModelBuilder();
				foreach (Vector v in tVectors){
					// builder.StartResidue();
					builder.AddAtom(new AtomInfo(Atom.NAME_C_ALPHA, Atom.ELEMENT_C, AtomTable.DEFAULT_ALT_LOC, false, v.X, v.Y, v.Z));
					builder.StartResidue();
				}
				new Protein(builder).SaveCif(Path.Combine(MainClass.Directory, "template-smooth.cif"));
				
				// new Protein (tVectors.Select ((v, i) => new Atom (i, Atom.NAME_C_ALPHA, ' ', "ALA", "X", i, ' ', v.X, v.Y, v.Z, 1, 1, "C", Atom.CHARGE_ZERO, false)))
				// 	.Save (Path.Combine (MainClass.Directory, "template-smooth.pdb"));

			}

			double scalingDistance = MainClass.STR_ALIGNMENT_SCALING_DISTANCE;
			double scalingFactor = 1.0 / (scalingDistance*scalingDistance);
			Func<int,int,double> scoreFunction1 = (i,j) => 1 / (1 + scalingFactor * (tVectors[i]-qVectors[j]).SqSize);
			Func<int,int,double> scoreFunction2 = (i, j) => LocalRMSDScore (tVectors, qVectors, i, j, 15, 7);
			Func<int,int,double> scoreFunction3 = (i, j) => scoreFunction1 (i, j) * scoreFunction2 (i, j);

			DateTime stamp = DateTime.Now;
			List<Tuple<int?,int?,double>> alignment = AlignPoints_DynProg (tVectors.Length, qVectors.Length, scoreFunction1);
			Lib.WriteLineDebug ("AlignPoints: {0}", DateTime.Now - stamp);

			if (Lib.DoWriteDebug) {
				StreamWriter w = new StreamWriter (Path.Combine (MainClass.Directory, "residue_alignment.tsv"));
				List<List<String>> alignmentTable = alignment.Select (t => new List<String> { 
					t.Item1 == null ? "-" : templates [(int)t.Item1].ChainId.ToString (), 
					t.Item1 == null ? "-" : templates [(int)t.Item1].SeqNumber.ToString (), 
					t.Item1 == null ? "-" : templates [(int)t.Item1].ShortName.ToString (), 
					t.Item3.ToString ("0.000"),
					t.Item2 == null ? "-" : candidates [(int)t.Item2].ShortName.ToString (),
					t.Item2 == null ? "-" : candidates [(int)t.Item2].SeqNumber.ToString (),
					t.Item2 == null ? "-" : candidates [(int)t.Item2].ChainId.ToString ()
				}).ToList ();
				w.WriteLine ("#chain1\tresi1\tresn1\tscore\tresn2\tresi2\tchain2");
				Lib.WriteCSV (w, alignmentTable, '\t');
				w.Close ();
			}

			return alignment;
		}

		private static double LocalRMSDScore(Vector[] template, Vector[] query, int i, int j, int window, double scalingR){
			double scalingFactor = 1 / scalingR;
			if (i < window || i + window >= template.Length || j < window || j + window >= query.Length) {
				return 0;
			} else {
				Matrix dump1;
				Matrix dump2;
				double rmsd;
				Matrix mobile = Matrix.FromRowVectors (query.GetRange (j - window, 2 * window + 1));
				Matrix target = Matrix.FromRowVectors (template.GetRange (i - window, 2 * window + 1));
				target.MeanCenterVertically (out dump1);
				LibAlgebra.Align (
					mobile,
					target,
					out dump1,
					out dump2,
					out rmsd
				);
				//LibAlgebra.FitRotation ();
				return 1 / (1 + scalingFactor*scalingFactor*rmsd*rmsd);
			}
		}

		private static Vector[] Smooth(Vector[] orig, int window){
			int n = orig.Length;
			Vector[] result = new Vector[n];
			double invCount = 1.0 / (2 * window + 1);
			for (int i = 0; i < window; i++) {
				result [i] = orig [i];
				result [n - i - 1] = orig [n - i - 1];
			}
			for (int i = window+1; i < n-window-1; i++) {
				result [i] = invCount * orig.GetRange (i - window, 2 * window + 1).Aggregate ((u, v) => u + v);
			}
			return result;
		}

		public static List<String> GetSequences(Chain ch, IEnumerable<SSE> sses){
			if (sses.Any (sse => sse.ChainID != ch.Id))
				throw new Exception ("Some SSEs belong to another chain");
			return ch.GetResidues (sses.Select (s => new Tuple<int,int> (s.Start, s.End)))
				.Select (s => String.Concat (s.Select (r => r.ShortName))).ToList ();
		}

		/*public static void ExtractSequences(String pdbFile, String annotationFile, String outputFile){
			Protein p = new Protein (pdbFile);
			List<SSE> sses = ReadAnnotationFile (annotationFile);
			IEnumerable<char> chains = sses.Where (sse=>!sse.IsNotFound ()).Select (sse=>sse.ChainID).Distinct ();

			List<Tuple<SSE,String>> ssesWithSequences;
			if (chains.Count () == 0) {
				ssesWithSequences = sses.Select (sse => new Tuple<SSE,String> (sse, "")).ToList ();
			} else if (chains.Count () == 1) {
				List<Tuple<int,int>> ranges = sses.Select (sse => sse.IsNotFound () ? new Tuple<int,int> (0, -1) : new Tuple<int,int> (sse.Start, sse.End)).ToList ();
				//ranges.Insert (0, new Tuple<int,int> (9, 12));
				List<List<Residue>> residues = p.GetChain (chains.First ()).GetResidues (ranges);
				ssesWithSequences = new List<Tuple<SSE,String>> ();
				for (int i = 0; i < sses.Count; i++) {
					String sequence = residues [i].Aggregate ("", (s, r) => s + r.ShortName);
					if (sses [i].Length () != sequence.Length) {
						Lib.WriteWarning ("This SSEs might contain missing residues: {0} (sequence: {1}).", sses [i], sequence);
					}
					ssesWithSequences.Add (new Tuple<SSE,String> (sses [i], sequence));
				}
			} else {
				throw new NotImplementedException (System.Reflection.MethodBase.GetCurrentMethod ().Name + ": This function is not implemented for SSEs from more than one chain.");
			}
			WriteAnnotationFileWithSequences (outputFile,ssesWithSequences,null);
		}*/
	
		public static List<List<Tuple<int,int>>> FindMaximalMatchings(int nG, List<Tuple<int,int>> edgesG,int nH, List<Tuple<int,int>> edgesH, Func<int,int,bool> subsequentInG, Func<int,int,bool> subsequentInH){
			
			//TODO put only pairs with low metric to modular product graph
			//TODO sort modular product graph vertices by metric
			//TODO use branch and bound in finding the best clique

			ISet<int> EG = new HashSet<int> (edgesG.Select (e => e.Item1 * nG + e.Item2)); //edge u->v is coded as u*nG+v
			ISet<int> EH = new HashSet<int> (edgesH.Select (e => e.Item1 * nG + e.Item2));
			ISet<int> subseqG = new HashSet<int> ();
			for (int u = 0; u < nG; u++)
				for (int v = 0; v < nG; v++)
					if (subsequentInG (u, v)) {
						subseqG.Add (u * nG + v);
						subseqG.Add (v * nG + u);
					}
			ISet<int> subseqH = new HashSet<int> ();
			for (int u = 0; u < nH; u++)
				for (int v = 0; v < nH; v++)
					if (subsequentInH (u, v)) {
						subseqG.Add (u * nH + v);
						subseqG.Add (v * nH + u);
					}
			
			int nM = nG * nH;
			bool[,] M = new bool[nM, nM];
			for (int u = 0; u < nG; u++) {
				for (int v = 0; v < nG; v++) {
					for (int u_ = 0; u_ < nH; u_++) {
						for (int v_ = 0; v_ < nH; v_++) {
							if (u!=v&&u_!=v_&&(u<v)==(u_<v_)&&
								(EG.Contains (u*nG+v)&&EH.Contains (u_*nH+v_)
								|| EG.Contains (v*nG+u)&&EH.Contains (v_*nH+u_)
								|| !EG.Contains (u*nG+v)&&!EG.Contains (v*nG+u)&&!EH.Contains (u_*nH+v_)&&!EH.Contains (v_*nH+u_)
								|| u==v && subseqH.Contains(u_*nH+v_)
								|| subseqG.Contains(u*nG+v_) && u_==v_
							)){
								M [u * nH + u_, v * nH + v_] = true; //tuple (u,u_) is coded as u*nH+u_
								M [v * nH + v_, u * nH + u_] = true;
							}
						}
					}
				}
			}
			List<List<int>> cliques = ListMaximalCliques (nM, M);
			return cliques.Select (c => c.Select (v => new Tuple<int,int> (v / nH, v % nH)).ToList ()).ToList ();
		}

		public static List<List<int>> ListMaximalCliques(int n, bool[,] neibMatrix){
			List<List<int>> cliques = new List<List<int>> ();
			ListMaximalCliques (n, neibMatrix, new List<int> (), 0, cliques);
			return cliques;
		}

		private static void ListMaximalCliques(int n, bool[,] neibMatrix, List<int> byNow, int firstAllowed, List<List<int>> output){
					bool largerFound = false;
			for (int i = firstAllowed; i < n; i++) {
				if (byNow.All (j => neibMatrix [i, j])) {
					largerFound = true;
					byNow.Add (i);
					ListMaximalCliques (n, neibMatrix, byNow, i + 1, output);
					byNow.RemoveAt (byNow.Count-1);
				}
			}
			if (!largerFound) {
				//this is a maximal clique
				if (output.All (c => !IsSubset (byNow, c)))
					output.Add (new List<int> (byNow));
			}
		}

		private static bool IsSubset(List<int> sub, List<int> super){
			return sub.Except (super).Count () == 0;
		}

		/** Finds common subgraph of graphs G and H with maximum sum of scores of matched vertices, keeping the order of vertices.
		 	Returns the list of matched vertex pairs.
			nG, nH - number of vertices of G, H
			edgesG, edgesH - adjacency matrices of G, H (0 = no edge, 1 = edge; if other values are used, they serve as edge labels - only edges with same label can be matched)
			score[i,j] - score for matching i-th vertex of G to j-th vertex of H*/
		public static List<Tuple<int,int>> MaxWeightOrderedMatching(int nG, int nH, int[,] edgesG, int[,] edgesH, double[,] score, bool[,] conflictsG, bool[,] conflictsH){
			if (edgesG.GetLength (0) != nG || edgesG.GetLength (1) != nG)
				throw new Exception ("Incorrect size of parameter edgesG.");
			if (edgesH.GetLength (0) != nH || edgesH.GetLength (1) != nH)
				throw new Exception ("Incorrect size of parameter edgesH.");
			if (score.GetLength (0) != nG || score.GetLength (1) != nH)
				throw new Exception ("Incorrect size of parameter score.");

			if (conflictsG==null) conflictsG = new bool[nG, nG];
			if (conflictsH==null) conflictsH = new bool[nH, nH];

			int nM = nG * nH;
			double[] w = Enumerable.Range (0, nM).Select (u => score [u / nH, u % nH]).ToArray ();
			Func<int,int,bool> edgesM = 
				(u, v) => 
				(u / nH == v / nH) == (u % nH == v % nH) //identity consistence
				&& (u / nH < v / nH) == (u % nH < v % nH) //order consistence
				&& edgesG [u / nH, v / nH] == edgesH [u % nH, v % nH] //connectivity consistence 
				&& !conflictsG [u/nH, v/nH] //vertices in G are not in conflict
				&& !conflictsH [u%nH, v%nH] //vertices in H are not in conflict
				&& w [u] > 0 //positive score for both matches
				&& w [v] > 0;
			int maxExpectedMatchingSize = Math.Min (nG, nH);

			List<int> maxWeightClique = MaxWeightClique (nM, edgesM, w, maxExpectedMatchingSize);
			return maxWeightClique.Select (u => new Tuple<int,int> (u / nH, u % nH)).ToList ();
		}

		public static List<Tuple<int,int>> MaxWeightMixedOrderedMatching(int nG, int nH, 
			int[] rank0G, int[] rank1G, int[] rank0H, int[] rank1H,
			double[,] score, bool softOrderConsistence /*, bool[,] conflictsG, bool[,] conflictsH*/)
		{

			if (rank0G.Length!= nG)
				throw new Exception ("Incorrect size of parameter rank0G.");
			if (rank1G.Length!= nG)
				throw new Exception ("Incorrect size of parameter rank1G.");
			if (rank0H.Length!= nH)
				throw new Exception ("Incorrect size of parameter rank0H.");
			if (rank1H.Length!= nH)
				throw new Exception ("Incorrect size of parameter rank1H.");
			if (score.GetLength (0) != nG || score.GetLength (1) != nH)
				throw new Exception ("Incorrect size of parameter score.");

			//if (conflictsG==null) conflictsG = new bool[nG, nG];
			//if (conflictsH==null) conflictsH = new bool[nH, nH];


			Func<int,int,int> ord; // order descriptor (allowed values: 0...k-1)
			bool[,] ordConsistencyMatrix; // order consistency matrix [k*k]
			if (softOrderConsistence){
				ord = (x, y) => Math.Min (Math.Max (-2, x - y), 2) + 2; // 5-state order descriptor (2 for equal, 1/3 for predecessor/successor, 0/4 for further less/greater)
				ordConsistencyMatrix = new int[,]{{1,1,0,0,0}, {1,1,1,0,0}, {0,1,1,1,0}, {0,0,1,1,1},{0,0,0,1,1}}.Select2D (x=>x!=0); // soft order consistency matrix
			} else {
				ord = (x, y) => Math.Sign (x - y) + 1; // 3-state order descriptor (1 for equal, 0 for less, 2 for greater)
				ordConsistencyMatrix = new int[,]{{1,0,0}, {0,1,0}, {0,0,1}}.Select2D (x=>x!=0); // hard order consistency matrix
			}
			int nM = nG * nH;
			double[] w = Enumerable.Range (0, nM).Select (u => score [u / nH, u % nH]).ToArray ();
			Func<int,int,bool> edgesM = 
				(u, v) => 
				(u / nH == v / nH) == (u % nH == v % nH) //identity consistence
				&& ordConsistencyMatrix[ord(rank0G[u/nH],rank0G[v/nH]),ord(rank0H[u%nH],rank0H[v%nH])] //order consistence of strands
				&& ordConsistencyMatrix[ord(rank0G[u/nH],rank1G[v/nH]),ord(rank0H[u%nH],rank1H[v%nH])]
				&& ordConsistencyMatrix[ord(rank1G[u/nH],rank0G[v/nH]),ord(rank1H[u%nH],rank0H[v%nH])]
				&& ordConsistencyMatrix[ord(rank1G[u/nH],rank1G[v/nH]),ord(rank1H[u%nH],rank1H[v%nH])]
				//&& !conflictsG [u/nH, v/nH] //vertices in G are not in conflict
				//&& !conflictsH [u%nH, v%nH] //vertices in H are not in conflict
				&& w [u] > 0 //positive score for both matches
				&& w [v] > 0;
			int maxExpectedMatchingSize = Math.Min (nG, nH);

			List<int> maxWeightClique = MaxWeightClique (nM, edgesM, w, maxExpectedMatchingSize);
			return maxWeightClique.Select (u => new Tuple<int,int> (u / nH, u % nH)).ToList ();
		}

		/** Returns inclusion-maximal clique with maximal total w. 
			Uses branch-and-bound algorithm.*/
		private static List<int> MaxWeightClique(int nM, Func<int,int,bool> edgesM, double[] w, int maxExpectedSize) {
			if (w.Length != nM)
				throw new Exception ("Incorrect size of parameter w.");

			int[] S = Enumerable.Range (0, nM).OrderBy (u => -w [u]).TakeWhile (u=>w[u]>0).ToArray (); //vertices ordered by descending weight, vertices with nonpositive weight are removed
			int nS = S.Count ();
			Lib.WriteLineDebug ("MaxWeightClique: {0} eligible vertices.", nS);
			List<int> K = new List<int>{ }; //best clique so far (indices of vertices in S)
			double wK = 0; //weight of best clique so far
			List<int> J = new List<int>{ }; //current clique (indices of vertices in S)
			double wJ = 0; //weight of current clique

			bool[,] Allowed = new bool[maxExpectedSize, nS]; //Allowed[J.Count,i] == S[i] is not in conflict with J; Allowed[J.Count,i] is valid only for future vertices i
			for (int i = 0; i < nS; i++) Allowed [0, i] = true;

			int p = 0; //index of current vertex
			bool forward = true; //direction of backtracking
			int iterationCounter = 0;

			while (p >= 0) {
				//Lib.WriteLineDebug ("{5} S[p]: {0}, Current: [{1}], ExpAdded: {2}, ExpScore: {3}, BestScore: {4}", p<nM?S [p]:-1, J.Select (u=>"("+u/12+","+u%12+")").EnumerateWithCommas (), maxExpAdded, maxExpScore, wK, forward?">":"<");
				if (forward) {
					if (p == nS || J.Count == maxExpectedSize) {
						//found new inclusion-maximal clique or all future vertices have nonpositive weights(useless to add them)
						//Lib.WriteLineDebug ("A");
						if (wJ > wK) {
							K = new List<int> (J);
							wK = wJ;
						}
						forward = false; p--;
					} else if (!Allowed [J.Count, p]) {
						//do not add vertex and go on forward
						//Lib.WriteLineDebug ("B");
						forward = true; p++;
					} else {
						int maxExpAdded = Math.Min (maxExpectedSize - J.Count, nS - p);
						double maxExpScore = wJ + Enumerable.Range (p, nS - p)
							.Where (i => Allowed [J.Count, i])
							.Take (maxExpAdded)
							.Select (i => w [S [i]])
							.Sum ();
						if (maxExpScore > wK) {
							//add new vertex and go on forward
							//Lib.WriteLineDebug ("C");
							J.Add (p);
							wJ = wJ + w [S [p]];
							if (J.Count < maxExpectedSize) {
								//prepare updated version of Allowed vector for the next iteration
								Array.Copy (Allowed, (J.Count - 1) * nS + p, Allowed, J.Count * nS + p, nS - p);
								for (int i = p + 1; i < nS; i++)
									if (!edgesM (S [p], S [i]))
										Allowed [J.Count, i] = false;
							}
							forward = true; p++;
						} else {
							//non-perspective branch
							//Lib.WriteLineDebug ("D");
							forward = false; p--;
						} 
					}
				} else { //backward
					if (J.Count > 0) {
						//remove trailing unused vertices and go forward again
						//Lib.WriteLineDebug ("E");
						p = J [J.Count - 1];
						J.RemoveAt (J.Count - 1);
						wJ = J.Sum (i => w [S [i]]);
						forward = true; p++;
					} else {
						//go backwards to the beginning and terminate algorithm
						//Lib.WriteLineDebug ("F");
						p = -1;
					}
				}
				iterationCounter++;
			}
			Lib.WriteLineDebug ("MaxWeightClique: {0} iterations.", iterationCounter);
			return K.Select (i => S [i]).ToList ();
		}

		/** Returns inclusion-maximal clique with maximal total w. 
			Uses branch-and-bound algorithm.*/
		private static List<int> MaxWeightClique_WithAllowedLists(int nM, Func<int,int,bool> edgesM, double[] w, int maxExpectedSize) {
			if (w.Length != nM)
				throw new Exception ("Incorrect size of parameter w.");

			int[] S = Enumerable.Range (0, nM).OrderBy (u => -w [u]).TakeWhile (u => w [u] > 0).ToArray (); //vertices ordered by descending weight, vertices with nonpositive weight are removed
			int nS = S.Count ();
			Console.WriteLine ("MaxWeightClique: {0} eligible vertices.", nS);
			if (nS == 0) {
				return new List<int> ();
			}

			int level = 0;
			int[] numAllowed = new int[maxExpectedSize]; // numAllowed[l] is the number of allowed vertices on level l
			numAllowed[0] = nS;
			int[,] allowed = new int[maxExpectedSize, nS]; // allowed[l, :] is the list of allowed vertices on level l (i.e. |current clique| = l) followed by a sentinel value nS
			for (int i = 0; i < nS; i++) {
				allowed [0, i] = S[i]; // on 0th level, all the vertices are allowed
			}
			int[] current = new int[maxExpectedSize+1];
			current [0] = -1;
			double[] weight = new double[maxExpectedSize+1]; // weight[l] is the weight on level l 
			weight[0] = 0;
			double[,] expWeight = new double[maxExpectedSize + 1, nS + 1];
			expWeight [0, 0] = 0;
			int[] bestClique = new int[maxExpectedSize]; // the best clique so far followed by a sentinel value nS
			int bestCliqueSize = 0;
			double bestWeight = 0;

			while (level >= 0) {
				//Console.WriteLine ($"Level {level}, current {current[level]}");
				if (level < maxExpectedSize && current[level] == -1){
					// I came from down
					for (int i = 0; i < numAllowed[level]; i++) {
						expWeight [level, i+1] = expWeight [level, i] + w [allowed [level, i]];
					}

					int maxExpAdded = Math.Min (numAllowed[level], maxExpectedSize - level);
					double maxExpWeight = weight [level];
					for (int i = 0; i < maxExpAdded; i++) {
						maxExpWeight += w [allowed [level, i]];
					}
					if (maxExpWeight <= bestWeight) {
						// go backwards
						level--;
						continue;
					}
				}

				// try extending
				current [level]++;
				if (level == maxExpectedSize || current[level] == numAllowed[level]){
					// used all allowed vertices on this level ==> not extensible
					if (weight[level] > bestWeight){
						// found new best clique
						bestCliqueSize = level;
						bestWeight = weight [level];
						for (int lev = 0; lev < level; lev++) {
							bestClique [lev] = allowed [lev, current [lev]];
						}
					}
					level--;
				} else {
					// extensible ==> add vertex current[level]
					int p = current [level];
					weight[level+1] = weight[level] + w[allowed[level, p]];
					if (level+1 < maxExpectedSize) {
						int nal = 0;
						for (int i = p+1; i < numAllowed[level]; i++) {
							if (edgesM(allowed[level, p], allowed[level, i])){
								allowed [level + 1, nal++] = allowed [level, i];
							}
						}
						numAllowed [level + 1] = nal;
					}
					current [level+1] = -1;
					level++;
				}
			}
			return bestClique.GetRange (0, bestCliqueSize).ToList ();
		}

	}
}


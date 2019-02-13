using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Cif.Components;

namespace protein
{
	public static class SecStrAssigners
	{
		public class SecStrAssignment{
			public List<SSE> SSEs{ get; set; }
			public List<Tuple<int,int,int>> Connectivity { get; set; } //each element is a triple <strand1,strand2,ladderType>, where ladderType = 1 for parallel, ladderType = -1 for antiparallel
			public List<Tuple<Residue,Residue>> HBonds{ get; set; }
			public List<Tuple<String,int,int>> MergeableSSEs{ get; set; } //each element is a triple <label,sse1,sse2> where sse1/sse2 are indices of the first/last SSE to be merged and label is the label of resulting merged SSE

			public SecStrAssignment(IEnumerable<SSE> sses){
				this.SSEs=sses.ToList ();
				this.Connectivity=new List<Tuple<int, int, int>>();
				this.HBonds=new List<Tuple<Residue, Residue>>();
				this.MergeableSSEs=new List<Tuple<String, int, int>>();
			}

			public static SecStrAssignment Combine(SecStrAssignment ass1, SecStrAssignment ass2){
				Lib.Shuffler shuffler;
				SecStrAssignment result = new SecStrAssignment (ass1.SSEs.ConcatAndGetShuffler (ass2.SSEs, out shuffler).ToList ());
				result.Connectivity = ass1.Connectivity.Concat (shuffler.UpdateIndices (ass2.Connectivity)).ToList ();
				result.HBonds = ass1.HBonds.Concat (ass2.HBonds).ToList ();
				return result;
			}
			public static SecStrAssignment Filter(SecStrAssignment ass, Func<SSE,bool> predicate){
				Lib.Shuffler shuffler;
				SecStrAssignment result = new SecStrAssignment (ass.SSEs.WhereAndGetShuffler (predicate, out shuffler));
				result.Connectivity = shuffler.UpdateIndices (ass.Connectivity).ToList ();
				result.HBonds = ass.HBonds;
				return result;
			}
			public static SecStrAssignment Order(SecStrAssignment ass){
				Lib.Shuffler shuffler;
				SecStrAssignment result = new SecStrAssignment (ass.SSEs.OrderAndGetShuffler (out shuffler));
				result.Connectivity = shuffler.UpdateIndices (ass.Connectivity).ToList ();
				result.HBonds = ass.HBonds;
				return result;
			}
		}

		public interface ISecStrAssigner {
			SecStrAssignment GetSecStrAssignment(); 
			String GetDescription();
		}

		public class SecStrAssignmentException : Exception{
			public SecStrAssignmentException(String message) : base(message) {}
			public SecStrAssignmentException(String message, Exception innerException) : base(message,innerException) {}
		}

		public class FileSecStrAssigner : ISecStrAssigner {
			public String SsaFile { get; private set; }
			public string[] ChainIDs { get; private set; }

			public FileSecStrAssigner(String ssaFile,  IEnumerable<string> chainIDs){
				this.SsaFile=ssaFile;
				this.ChainIDs=chainIDs.ToArray ();
			}

			public SecStrAssignment GetSecStrAssignment(){
				try {
					return new SecStrAssignment(LibAnnotation.ReadAnnotationFile (SsaFile).Where (x => ChainIDs.Contains (x.ChainID)).ToList ());
				} catch (IOException e) {
					Lib.WriteError ("Could not open \"" + SsaFile + "\".");
					throw new SecStrAssignmentException (MethodBase.GetCurrentMethod ().DeclaringType.Name + " failed.", e);
				}
			}

			public String GetDescription(){
				return "reading from file " + Path.GetFileName (SsaFile) + " (chainIDs: " + ChainIDs.EnumerateWithCommas () + ")";
			}
		}

		public class FileSecStrAssigner_Json : ISecStrAssigner {
			public String SsaFile { get; private set; }
			public String EntryName{ get; private set; }
			public string[] ChainIDs { get; private set; }

			/* If entryName==null, it will take the first entry in the file. */
			public FileSecStrAssigner_Json(String ssaFile,String entryName, IEnumerable<string> chainIDs){
				this.SsaFile=ssaFile;
				this.EntryName=entryName;
				this.ChainIDs=chainIDs.ToArray ();
			}

			public SecStrAssignment GetSecStrAssignment(){
				try {
					String dump;
					List<Tuple<int,int,int>> connectivity;
					List<Tuple<String,int,int>> mergeable;
					Lib.Shuffler shuffler;
					SecStrAssignment result = new SecStrAssignment (LibAnnotation.ReadAnnotationFile_Json (SsaFile,EntryName,out dump,out connectivity,out mergeable, true).WhereAndGetShuffler (x => ChainIDs.Contains (x.ChainID), out shuffler).ToList ());
					connectivity=shuffler.UpdateIndices (connectivity).ToList ();
					result.Connectivity=connectivity;
					result.MergeableSSEs=mergeable;
					return result;
				} catch (FormatException e){
					Lib.WriteError ("Invalid format of the template annotation file:\n    " + e.Message );
					throw new SecStrAssignmentException (MethodBase.GetCurrentMethod ().DeclaringType.Name + " failed.", e);
				} catch (IOException e) {
					Lib.WriteError ("Could not open \"" + SsaFile + "\".");
					throw new SecStrAssignmentException (MethodBase.GetCurrentMethod ().DeclaringType.Name + " failed.", e);
				}
			}

			public String GetDescription(){
				return "reading from JSON file " + Path.GetFileName (SsaFile) + " (chainIDs: " + ChainIDs.EnumerateWithCommas () + ")";
			}
		}

		public class DsspSecStrAssigner : ISecStrAssigner {
			public String DSSPExecutable { get; private set; }
			public String PDBFile { get; private set; }
			public String DSSPFile { get; private set; }
			public string[] ChainIDs { get; private set; }
			public char[] AcceptedSSETypes { get; private set; }

			public DsspSecStrAssigner(String dsspExecutable, String PDBFile, String DSSPFile, IEnumerable<string> chainIDs, char[] acceptedSSETypes){
				this.DSSPExecutable = dsspExecutable;
				this.PDBFile = PDBFile;
				this.DSSPFile = DSSPFile;
				this.ChainIDs=chainIDs.ToArray ();
				this.AcceptedSSETypes = acceptedSSETypes;
			}

			public SecStrAssignment GetSecStrAssignment(){
				Lib.WriteInColor (ConsoleColor.Yellow, "Running DSSP:\n");
				if (!Lib.RunDSSP (DSSPExecutable, PDBFile, DSSPFile))
					throw new  SecStrAssignmentException (MethodBase.GetCurrentMethod ().DeclaringType.Name + " failed.");
				
				try { 
					return new SecStrAssignment(LibAnnotation.ReadSSEsFromDSSP (DSSPFile, AcceptedSSETypes)
						.Where (x => ChainIDs.Contains (x.ChainID))
						.Select ((sse, i) => sse.RelabeledCopy (sse.Type + i.ToString ()))
						.ToList ()); 
				}
				catch (Exception e) { throw new SecStrAssignmentException (MethodBase.GetCurrentMethod ().DeclaringType.Name + " failed.", e); }
			}

			public String GetDescription(){
				return "DSSP method for " + Path.GetFileName (PDBFile) + " (accepted types: " + Lib.EnumerateWithCommas (AcceptedSSETypes) + ")";
			}
		}

		public class SingleChainGeomSecStrAssigner : ISecStrAssigner {
			public Chain Chain { get; private set; }
			public double RmsdLimit { get; private set; }
			public bool HelicesAllowed { get; private set; }
			public bool SheetsAllowed { get; private set; }

			public SingleChainGeomSecStrAssigner(Chain chain, double rmsdLimit){
				this.Chain=chain;
				this.RmsdLimit=rmsdLimit;
				this.HelicesAllowed=true;
				this.SheetsAllowed=true;
			}

			public SingleChainGeomSecStrAssigner(Chain chain, double rmsdLimit, char[] acceptedSSETypes){
				this.Chain=chain;
				this.RmsdLimit=rmsdLimit;
				this.HelicesAllowed= acceptedSSETypes.Contains (SSE.MIXED_HELIX_TYPE) || acceptedSSETypes.Contains (SSE.HELIX_H_TYPE);
				this.SheetsAllowed=acceptedSSETypes.Contains (SSE.SHEET_TYPE);
			}

			public SecStrAssignment GetSecStrAssignment(){
				int firstAssignedInQuad = 1;
				int lastAssignedInQuad = 2;
				int minimumSSELength = 3;
				List<SSE> result = new List<SSE> ();

				List<Residue> residues = Chain.GetResidues ().Where (r => r.HasCAlpha ()).ToList ();
				bool makingHelix = false;
				bool makingSheet = false;
				char helixType = SSE.MIXED_HELIX_TYPE;
				char sheetType = SSE.SHEET_TYPE;
				int start = 0;
				int end = 0;
				double rmsd = 0;
				double maxRmsd = 0;
				for (int i = 0; i < residues.Count-3; i++) {
					List<Residue> quad = residues.GetRange (i, 4);
					//Console.WriteLine ("at {0}",residues[i].ResSeq);
					if (makingHelix) {
						if (LibAnnotation.CheckGeometryOf1Unit (quad, helixType, RmsdLimit, out rmsd) == LibAnnotation.GeometryCheckResult.OK) {
							//continue helix
							end = i + lastAssignedInQuad;
							maxRmsd = Math.Max (maxRmsd, rmsd);
						} else {
							//finish helix
							if (residues [end].SeqNumber - residues [start].SeqNumber + 1 >= minimumSSELength)
								result.Add (new SSE ("_", Chain.Id, residues [start].SeqNumber, residues [end].SeqNumber, helixType,null));
							makingHelix = false;
						}
					} else if (makingSheet) {
						if (LibAnnotation.CheckGeometryOf1Unit (quad, sheetType, RmsdLimit, out rmsd) == LibAnnotation.GeometryCheckResult.OK) {
							//continue sheet
							end = i + lastAssignedInQuad;
							maxRmsd = Math.Max (maxRmsd, rmsd);
						} else {
							//finish sheet
							//Console.WriteLine ("Ending sheet at {0}", residues [i + lastAssignedInQuad].ResSeq);
							if (residues [end].SeqNumber - residues [start].SeqNumber + 1 >= minimumSSELength)
								result.Add (new SSE ("_", Chain.Id, residues [start].SeqNumber, residues [end].SeqNumber, sheetType,null));
							makingSheet = false;
						}
					} 
					if (!makingHelix && !makingSheet) {
						if (HelicesAllowed && LibAnnotation.CheckGeometryOf1Unit (quad, helixType, RmsdLimit, out rmsd) == LibAnnotation.GeometryCheckResult.OK) {
							//start helix
							//Console.WriteLine ("Starting helix at {0}", residues [i + firstAssignedInQuad].ResSeq);
							makingHelix = true;
							start = i + firstAssignedInQuad;
							end = i + lastAssignedInQuad;
							maxRmsd = rmsd;
						} else if (SheetsAllowed && LibAnnotation.CheckGeometryOf1Unit (quad, sheetType, RmsdLimit, out rmsd) == LibAnnotation.GeometryCheckResult.OK) {
							//start sheet
							//Console.WriteLine ("Starting sheet at {0}", residues [i + firstAssignedInQuad].ResSeq);
							makingSheet = true;
							start = i + firstAssignedInQuad;
							end = i + lastAssignedInQuad;
							maxRmsd = rmsd;
						}
					}
				}
				if (makingHelix) {
					//finish helix
					if (residues [end].SeqNumber - residues [start].SeqNumber + 1 >= minimumSSELength)
						result.Add (new SSE ("_", Chain.Id, residues [start].SeqNumber, residues [end].SeqNumber, helixType,null));
					makingHelix = false;
				}
				if (makingSheet){
					//finish sheet
					if (residues [end].SeqNumber - residues [start].SeqNumber + 1 >= minimumSSELength)
						result.Add (new SSE ("_", Chain.Id, residues [start].SeqNumber, residues [end].SeqNumber, sheetType,null));
					makingSheet = false;
				}

				if (SheetsAllowed) {
					// Deleting sheets that have no hydrogen bonds
					List<List<Residue>> sheets = Chain.GetResidues (result.Where (sse => sse.Type == sheetType).Select (sse => new Tuple<int,int> (sse.Start, sse.End)));
					IEnumerable<Atom> atomsN = sheets.SelectMany (sse => sse.SelectMany (r => r.GetNAmides()));
					IEnumerable<Atom> atomsO = sheets.SelectMany (sse => sse.SelectMany (r => r.GetOCarbs()));
					double maxDistanceForHydrogenBond = 3.5; //TODO Get a rational value for this number (distance N-O). This value is just a guess.
					int minResiduesBetweenHydrogenBond = 2; // number of residues that must be between 2 residues that are forming a beta-sheet stabilizing hydrogen bond

					List<int> deletedSheetsByStart = new List<int> ();
					for (int i = sheets.Count - 1; i >= 0; i--) {
						int firstResSeq = sheets [i] [0].SeqNumber;
						int lastResSeq = sheets [i] [sheets [i].Count - 1].SeqNumber;

						// Version without shortening sheets
						/*var atomsNHere = sheets [i].SelectMany (r => r.GetAtoms ().Where (a => a.Name == " N  "));
				var atomsOHere = sheets [i].SelectMany (r => r.GetAtoms ().Where (a => a.Name == " O  "));
				if (!atomsNHere.Any (a1 => atomsO.Any (a2 => !(a2.ResSeq >= firstResSeq && a2.ResSeq <= lastResSeq) && (Math.Abs (a2.ResSeq-a1.ResSeq)>minResiduesBetweenHydrogenBond) && LibProtein.Distance (a1, a2) <= maxDistanceForHydrogenBond))
					&& !atomsOHere.Any (a1 => atomsN.Any (a2 => !(a2.ResSeq >= firstResSeq && a2.ResSeq <= lastResSeq) && (Math.Abs (a2.ResSeq-a1.ResSeq)>minResiduesBetweenHydrogenBond) && LibProtein.Distance (a1, a2) <= maxDistanceForHydrogenBond))) {
					deletedSheetsByStart.Add (firstResSeq);
					Console.WriteLine ("Deleting sheet at {0}.",firstResSeq);
				}else
					Console.WriteLine ("Keeping sheet at {0}.",firstResSeq);
				result.RemoveAll (sse => sse.Type == sheetType && deletedSheetsByStart.Contains (sse.Start));*/

						// Version with shortening sheets
						List<int> hbN = sheets [i].Where (r => r.GetNAmides().Any (a1 => atomsO.Any (a2 => 
							!(a2.ResidueSeqNumber >= firstResSeq && a2.ResidueSeqNumber <= lastResSeq)
							&& (Math.Abs (a2.ResidueSeqNumber - a1.ResidueSeqNumber) > minResiduesBetweenHydrogenBond)
							&& LibProtein.Distance (a1, a2) <= maxDistanceForHydrogenBond))).Select (r => r.SeqNumber).ToList ();
						List<int> hbO = sheets [i].Where (r => r.GetOCarbs().Any (a1 => atomsN.Any (a2 => 
							!(a2.ResidueSeqNumber >= firstResSeq && a2.ResidueSeqNumber <= lastResSeq)
							&& (Math.Abs (a2.ResidueSeqNumber - a1.ResidueSeqNumber) > minResiduesBetweenHydrogenBond)
							&& LibProtein.Distance (a1, a2) <= maxDistanceForHydrogenBond))).Select (r => r.SeqNumber).ToList ();
						List<int> hb = hbN.Union (hbO).ToList ();
						if (hb.Count == 0) {
							result.RemoveAll (sse => sse.Start == firstResSeq);
						} else {
							int index = result.FindIndex (sse => sse.Start == firstResSeq);
							SSE old = result [index];
							//result [index] = new SSE (old.Label, old.ChainID, Math.Max (old.Start, hb.Min ()-1), Math.Min (old.End, hb.Max ()+1), result [index].Type);
							result [index] = new SSE (old.Label, old.ChainID, hb.Min (), hb.Max (), result [index].Type,null);
						}
					}
				}

				return new SecStrAssignment(result.Select ((sse, i) => sse.RelabeledCopy (sse.Type + i.ToString ())).ToList ());
			}

			public String GetDescription(){
				return "geometrical method with RMSD limit " + RmsdLimit + " A (accepted types: " + (HelicesAllowed?SheetsAllowed?"H, E":"H":SheetsAllowed?"E":"" ) + ")";
			}
		}

		public class GeomSecStrAssigner : ISecStrAssigner {
			public ISecStrAssigner[] Assigners { get; private set; }

			public GeomSecStrAssigner(IEnumerable<Chain> chains, double rmsdLimit){
				Assigners=chains.Select (c=>new SingleChainGeomSecStrAssigner (c,rmsdLimit)).ToArray ();
			}

			public GeomSecStrAssigner(IEnumerable<Chain> chains, double rmsdLimit, char[] acceptedSSETypes){
				Assigners=chains.Select (c=>new SingleChainGeomSecStrAssigner (c,rmsdLimit,acceptedSSETypes)).ToArray ();
			}

			public SecStrAssignment GetSecStrAssignment(){
				SecStrAssignment result = new SecStrAssignment (new List<SSE>());
				foreach(ISecStrAssigner assigner in Assigners){
					result = SecStrAssignment.Combine (result, assigner.GetSecStrAssignment ());
				}
				return SecStrAssignment.Order (result);
			}

			public String GetDescription(){
				return Assigners.Count () > 0 ? Assigners [0].GetDescription () : "";
			}
		}

		public class GeomDsspSecStrAssigner : ISecStrAssigner {
			public ISecStrAssigner HelixAssigner { get; private set; }
			public ISecStrAssigner SheetAssigner { get; private set; }

			public GeomDsspSecStrAssigner(ISecStrAssigner helixAssigner, ISecStrAssigner sheetAssigner){
				this.HelixAssigner=helixAssigner;
				this.SheetAssigner=sheetAssigner;
			}

			public GeomDsspSecStrAssigner(IEnumerable<Chain> chains, double rmsdLimit, String dsspExecutable, String PDBFile, String DSSPFile, char[] acceptedSSETypes){
				this.HelixAssigner = new GeomSecStrAssigner(chains, rmsdLimit, acceptedSSETypes.Intersect (SSE.ALL_HELIX_TYPES).ToArray ());
				this.SheetAssigner = new DsspSecStrAssigner(dsspExecutable, PDBFile, DSSPFile, chains.Select (c=>c.Id), acceptedSSETypes.Intersect (SSE.ALL_SHEET_TYPES).ToArray ());
			}

			public SecStrAssignment GetSecStrAssignment(){
				SecStrAssignment result = SecStrAssignment.Order (SecStrAssignment.Combine (SheetAssigner.GetSecStrAssignment (), HelixAssigner.GetSecStrAssignment ()));
				//reporting overlapping SSEs
				for (int i = 0; i < result.SSEs.Count - 1; i++) {
					if (result.SSEs [i].End >= result.SSEs [i + 1].Start) {
						Lib.WriteLineDebug ("Overlapping SSEs: {0} and {1}", result.SSEs [i], result.SSEs [i + 1]);
					}
				}
				// truncating overlapping SSEs
				/*if (result.Select (sse=>sse.ChainID).Distinct ().Count () > 1) 
					throw new NotImplementedException ("This method is not implemented for SSEs coming from multiple chains.");
				for (int i = 0; i < result.Count-1; i++) {
					if (result [i].End >= result [i + 1].Start){
						//TODO This solution might not be perfect.
						result[i].End = result[i+1].Start-1;
						Lib.WriteLineDebug ("Some SSEs were overlapping and therefore they have been truncated.");
					}
				}*/
				return result;
			}

			public String GetDescription(){
				return "mixed method (for helices: " + HelixAssigner.GetDescription () + ", for sheets: " + SheetAssigner.GetDescription ();
			}
		}

		public class GeomHbondSecStrAssigner : ISecStrAssigner {
			public ISecStrAssigner GeomAssigner { get; private set; }
			public HBondSecStrAssigner HBondAssigner { get; private set; }
			public const int MIN_HELIX_SUBHELIX_OVERLAP = 1;

			public GeomHbondSecStrAssigner(Protein p, double rmsdLimit, double hBondEnergyCutoff){
				this.GeomAssigner = new GeomSecStrAssigner(p.GetChains (), rmsdLimit, SSE.ALL_HELIX_TYPES);
				this.HBondAssigner = new HBondSecStrAssigner(p, hBondEnergyCutoff);
			}

			public SecStrAssignment GetSecStrAssignment(){
				HBondAssigner.DetectSheets = true;
				HBondAssigner.DetectHelices = false;
				SecStrAssignment sheetAssignment = HBondAssigner.GetSecStrAssignment ();
				HBondAssigner.DetectSheets = false;
				HBondAssigner.DetectHelices = true;
				List<SSE> subhelices = HBondAssigner.GetSecStrAssignment ().SSEs.OrderBy (sse => sse).ToList ();
				HBondAssigner.DetectSheets = true;
				List<SSE> bigHelices = GeomAssigner.GetSecStrAssignment ().SSEs.OrderBy (sse => sse).ToList ();

				Func<SSE,SSE,int,bool> doOverlap = (a, b, mo) => a.ChainID == b.ChainID && a.End >= b.Start + mo && b.End >= a.Start + mo;

				List<SSE> processedHelices = new List<SSE> ();
				foreach (var helix in bigHelices) {
					List<SSE> subs = subhelices.Where (sub => doOverlap (helix, sub, MIN_HELIX_SUBHELIX_OVERLAP)).ToList ();
					if (subs.Count > 0) {
						foreach (var sub in subs) {
							helix.AddNestedSSE (sub);
						}
						var subTypes = subs.Select (sub => sub.Type).Distinct ();
						helix.Type = subTypes.Count () == 1 ? subTypes.First () : SSE.MIXED_HELIX_TYPE;
						processedHelices.Add (helix);
					} else {
						// ignore helix that contains no subhelices (no H-bonds)
					}
				}
				var result = SecStrAssignment.Order(SecStrAssignment.Combine (sheetAssignment, new SecStrAssignment (processedHelices)));

				//reporting overlapping SSEs
				var sses = result.SSEs.OrderBy (x=>x).ToArray ();
				for (int i = 0; i < sses.Length - 1; i++) {
					if (doOverlap(sses[i], sses[i+1], 0))
						Lib.WriteLineDebug ("Overlapping SSEs: {0} and {1}", sses [i], sses [i + 1]);
				}

				return result;
			}

			public String GetDescription(){
				return "mixed method (for helices: " + GeomAssigner.GetDescription () + " with further division using " + GeomAssigner.GetDescription () + ", for sheets: " + HBondAssigner.GetDescription ();
			}
		}

		public class JoiningSecStrAssigner : ISecStrAssigner {
			public ISecStrAssigner InnerAssigner { get; private set; }
			public Protein Protein { get; private set; }
			public double JoiningRmsdLimit { get; private set; }
			public Func<char,char,char?> JoiningTypeCombining { get; private set; }

			public JoiningSecStrAssigner(ISecStrAssigner innerAssigner, Protein protein, double joiningRmsdLimit, Func<char,char,char?> joiningTypeCombining){
				this.InnerAssigner = innerAssigner;
				this.Protein = protein;
				this.JoiningRmsdLimit = joiningRmsdLimit;
				this.JoiningTypeCombining = joiningTypeCombining;
			}
			

			public SecStrAssignment GetSecStrAssignment(){
				//TODO connectivity
				throw new NotImplementedException ();
				SecStrAssignment unjoined = InnerAssigner.GetSecStrAssignment ();
				return new SecStrAssignment (LibAnnotation.JoinSSEs_GeomVersion (Protein, unjoined.SSEs, JoiningRmsdLimit, JoiningTypeCombining));
			}

			public String GetDescription(){
				return InnerAssigner.GetDescription () + " and then joining using geometric criterion with RMSD limit " + JoiningRmsdLimit.ToString ("0.000") + " A";
			}
		}

		public class OutputtingSecStrAssigner : ISecStrAssigner {
			public ISecStrAssigner InnerAssigner { get; private set; }
			public String OutputFile { get; private set; }

			public OutputtingSecStrAssigner(ISecStrAssigner innerAssigner, String outputFile){
				this.InnerAssigner = innerAssigner;
				this.OutputFile = outputFile;
			}

			public SecStrAssignment GetSecStrAssignment(){
				SecStrAssignment result = InnerAssigner.GetSecStrAssignment ();
				LibAnnotation.WriteAnnotationFile (OutputFile, result.SSEs, "Helix info obtained by " + InnerAssigner.GetDescription() + ".");
				//TODO write out connectivity and hbonds
				return result;
			}

			public String GetDescription(){
				return InnerAssigner.GetDescription ();
			}
		}

		public class OutputtingSecStrAssigner_Json : ISecStrAssigner {
			public ISecStrAssigner InnerAssigner { get; private set; }
			public String EntryName{ get; private set; }
			public String OutputFile { get; private set; }

			public OutputtingSecStrAssigner_Json(ISecStrAssigner innerAssigner, String outputFile, String entryName){
				if (entryName==null) throw new ArgumentNullException ();
				this.InnerAssigner = innerAssigner;
				this.EntryName=entryName;
				this.OutputFile = outputFile;
			}

			public SecStrAssignment GetSecStrAssignment(){
				SecStrAssignment result = InnerAssigner.GetSecStrAssignment ();
				LibAnnotation.WriteAnnotationFile_Json (OutputFile, 
					EntryName, 
					result.SSEs,
					null, 
					result.Connectivity,
					result.HBonds,
					"Helix info obtained by " + InnerAssigner.GetDescription() + "."
				);
				//TODO write out connectivity and hbonds
				return result;
			}

			public String GetDescription(){
				return InnerAssigner.GetDescription ();
			}
		}

		public class RelabellingSecStrAssigner : ISecStrAssigner {
			/*Returns the same SSEs as its InnerAssigner, but relabels them with this scheme:
				<Prefix><SSE_type><sequentially_assigned_identifier>, e.g. 1tqn_H10 */
			public ISecStrAssigner InnerAssigner { get; private set; }
			public String Prefix { get; private set; }
			public bool SetAllLabelsToNull{ get; private set; }

			public RelabellingSecStrAssigner(ISecStrAssigner innerAssigner, String prefix, bool setAllLabelToNull){
				this.InnerAssigner = innerAssigner;
				this.Prefix = prefix;
				this.SetAllLabelsToNull=setAllLabelToNull;
			}

			public SecStrAssignment GetSecStrAssignment(){
				SecStrAssignment result = InnerAssigner.GetSecStrAssignment ();
				result.SSEs=result.SSEs
					.Select ((sse, i) => sse.RelabeledCopy (SetAllLabelsToNull?null: Prefix + sse.Type + i.ToString ()))
					.ToList ();
				return result;
			}

			public String GetDescription(){
				return InnerAssigner.GetDescription ();
			}
		}

		public class FilteringSecStrAssigner : ISecStrAssigner {
			/*Returns the same SSEs as its InnerAssigner, but relabels them with this scheme:
				<Prefix><SSE_type><sequentially_assigned_identifier>, e.g. 1tqn_H10 */
			public ISecStrAssigner InnerAssigner { get; private set; }
			public char[] AcceptedSSETypes { get; private set; }
			public string[] ChainIDs{ get; private set; }

			public FilteringSecStrAssigner(ISecStrAssigner innerAssigner, IEnumerable<char> acceptedTypes, IEnumerable<string> chainIDs){
				this.InnerAssigner = innerAssigner;
				this.AcceptedSSETypes=acceptedTypes.ToArray ();
				this.ChainIDs=chainIDs.ToArray ();
			}

			public SecStrAssignment GetSecStrAssignment(){
				/*Lib.Shuffler shuffler;
				List<SSE> result = InnerAssigner
					.GetSecStrAssignment (out connectivity)
					.WhereAndGetShuffler (sse => AcceptedSSETypes.Contains (sse.Type) && ChainIDs.Contains (sse.ChainID), out shuffler)
					.ToList ();
				connectivity=shuffler.UpdateIndices (connectivity).ToList ();*/
				return SecStrAssignment.Filter (InnerAssigner.GetSecStrAssignment (), sse => AcceptedSSETypes.Contains (sse.Type) && ChainIDs.Contains (sse.ChainID));
			}

			public String GetDescription(){
				return InnerAssigner.GetDescription ();
			}
		}

		private interface IHBondFinder{
			bool IsHBond (int donor, int acceptor);
			List<int> FindHAcceptors(int donor);
			List<int> FindHDonors(int acceptor);
		}

		private class SimpleHBondFinder : IHBondFinder
		{
			protected Residue[] residues;
			private Vector[] vecH;
			private Vector[] vecN;
			protected Vector[] vecCA;
			private Vector[] vecC;
			private Vector[] vecO;
			private bool[] canBeDonor;
			private bool[] canBeAcceptor;
			private double energyCutoff;

			public SimpleHBondFinder(IEnumerable<Residue> residues, double dsspEnergyCutoff){
				this.residues = residues.ToArray();
				energyCutoff=dsspEnergyCutoff;

				canBeDonor=new bool[this.residues.Length];
				canBeAcceptor=new bool[this.residues.Length];
				vecCA=new Vector[this.residues.Length];
				vecH=new Vector[this.residues.Length];
				vecN=new Vector[this.residues.Length];
				vecC=new Vector[this.residues.Length];
				vecO=new Vector[this.residues.Length];
				for (int i = 0; i < this.residues.Length; i++) {
					Residue r = this.residues[i];
					Atom? cAlpha = r.GetCAlpha();
					if (cAlpha == null){
						throw new ArgumentException ($"Residue {r} has no C-alpha atom.");
					}
					vecCA[i]= ((Atom) r.GetCAlpha()).Position ();
					Atom? hAmide = r.GetHAmide();
					Atom? nAmide = r.GetNAmide();
					Atom? cCarb = r.GetCCarb();		
					Atom? oCarb = r.GetOCarb();				
					if (hAmide != null && nAmide != null){
						canBeDonor[i]=true;
						vecH[i]=((Atom) hAmide).Position ();
						vecN[i]=((Atom) nAmide).Position ();
					}
					if (cCarb != null && oCarb != null){
						canBeAcceptor[i]=true;
						vecC[i]=((Atom) cCarb).Position ();
						vecO[i]=((Atom) oCarb).Position ();
					}
				}
			}

			public double DsspEnergy(int donor, int acceptor){
				if (canBeDonor [donor] && canBeAcceptor [acceptor]) {
					Vector n = vecN [donor];
					Vector h = vecH [donor];
					Vector c = vecC [acceptor];
					Vector o = vecO [acceptor];
					return 0.084 * 332 * (1 / (o - n).Size + 1 / (c - h).Size - 1 / (o - h).Size - 1 / (c - n).Size);
				} else {
					throw new InvalidOperationException ("Some residues miss atoms which are important for calculation of DSSP energy.");
				}
			}

			public bool IsHBond(int donor, int acceptor){
				if (donor < 0 || acceptor < 0)
					return false;
				if (residues [donor].ChainId == residues [acceptor].ChainId && Math.Abs (residues [donor].SeqNumber - residues [acceptor].SeqNumber) <= 1) {
					return false;
				}
				try {
					return DsspEnergy (donor,acceptor) <= energyCutoff;
				} catch (InvalidOperationException) {
					//some residues do not contain needed atoms
					return false;
				}
			}

			public virtual List<int> FindHAcceptors(int res){
				return Enumerable.Range (0, residues.Length).Where (s => IsHBond (res, s)).ToList ();
			}

			public virtual List<int> FindHDonors(int res){
				return Enumerable.Range (0, residues.Length).Where (s => IsHBond (s, res)).ToList ();
			}
		}

		private class BoxingHBondFinder:SimpleHBondFinder{
			private const double BOX_SIZE=10;
			private List<int>[,,] boxes;
			private List<int>[,,] neighborBoxes;
			private double[] corner;
			private int[] nBoxes;
			public BoxingHBondFinder(IEnumerable<Residue> residues_, double dsspEnergyCutoff)
				: base(residues_,dsspEnergyCutoff){
				corner=new double[]{vecCA.Select (v=>v.X).Min (),vecCA.Select (v=>v.Y).Min (),vecCA.Select (v=>v.Z).Min ()};
				double[] otherCorner=new double[]{vecCA.Select (v=>v.X).Max (),vecCA.Select (v=>v.Y).Max (),vecCA.Select (v=>v.Z).Max ()};
				nBoxes=new int[3];
				for (int i=0;i<3;i++) nBoxes[i]=(int)Math.Floor((otherCorner[i]-corner[i])/BOX_SIZE)+1;
				boxes=new List<int>[nBoxes[0],nBoxes[1],nBoxes[2]];
				for (int i = 0; i < nBoxes[0]; i++) {
					for (int j = 0; j < nBoxes[1]; j++) {
						for (int k = 0; k < nBoxes[2]; k++) {
							boxes[i,j,k]=new List<int>();
						}
					}
				}
				for (int i=0; i<vecCA.Length;i++) {
					GetBox (i).Add (i);
				}
				neighborBoxes=new List<int>[nBoxes[0],nBoxes[1],nBoxes[2]];
				for (int i = 0; i < nBoxes[0]; i++) {
					for (int j = 0; j < nBoxes[1]; j++) {
						for (int k = 0; k < nBoxes[2]; k++) {
							neighborBoxes[i,j,k]=MergeNeighborBoxes (new int[]{i,j,k});
						}
					}
				}
			}
			private int[] GetBoxIndices(int res){
				int[] indices = new int[3];
				double[] pos = vecCA[res].AsArray ();
				for (int i = 0; i < 3; i++)
					indices [i] = (int)Math.Floor ((pos [i] - corner [i]) / BOX_SIZE);
				return indices;
			}
			private List<int> GetBox(int res){
				int[] indices = GetBoxIndices (res);
				return boxes [indices [0], indices [1], indices [2]];
			}
			private List<int> GetNeighborBox(int res){
				int[] indices = GetBoxIndices (res);
				return neighborBoxes [indices [0], indices [1], indices [2]];
			}
			private List<int> MergeNeighborBoxes(int[] indices){
				List<List<int>> result = new List<List<int>> ();
				for (int i = Math.Max (0,indices[0]-1); i <= Math.Min (nBoxes[0]-1,indices[0]+1); i++) {
					for (int j = Math.Max (0, indices [1] - 1); j <= Math.Min (nBoxes [1]-1, indices [1] + 1); j++) {
						for (int k = Math.Max (0, indices [2] - 1); k <= Math.Min (nBoxes [2]-1, indices [2] + 1); k++) {
							result.Add (boxes [i, j, k]);
						}
					}
				}
				return result.SelectMany (b => b).ToList ();
			}

			public override List<int> FindHAcceptors(int res){
				return GetNeighborBox (res).Where (s => IsHBond (res, s)).ToList ();
			}

			/*public virtual List<int> FindHDonors(Residue r){
				return Lib.IndicesWhere(residues, s => IsHBond(s,r));
			}*/

			public override List<int> FindHDonors(int res){
				return GetNeighborBox (res).Where (s => IsHBond (s, res)).ToList ();
			}

		}

		private class BetaLadder{
			public enum LadderType { Parallel, Antiparallel }
			public LadderType Type { get; private set; }

			public enum HBondDirection { From0To1, From1To0 }
			private HBondDirection InvertedDirection(HBondDirection d){
				return d == HBondDirection.From0To1 ? HBondDirection.From1To0 : HBondDirection.From0To1;
			}
			public HBondDirection FirstHBondDirection { get; set; }
			public HBondDirection LastHBondDirection { get; set; }

			public int Start0 { get; set; }
			public int End0 { get; set; }

			public int Start1 { get; set; }
			public int End1 { get; set; }

			// Z is defined as follows: Z=3*resi for backbone NH group, Z=3*resi+2 for backbone CO group
			public int ZStart0{ get { return 3 * Start0 + (FirstHBondDirection == HBondDirection.From1To0 ? 2 : 0); } }
			public int ZEnd0{ get { return 3 * End0 + (LastHBondDirection == HBondDirection.From1To0 ? 2 : 0); } }
			public int ZStart1 {
				get { 
					if (Type == LadderType.Parallel)
						return 3 * Start1 + (FirstHBondDirection == HBondDirection.From0To1 ? 2 : 0);
					else
						return 3 * Start1 + (LastHBondDirection == HBondDirection.From0To1 ? 2 : 0); 
				}
			}
			public int ZEnd1 {
				get { 
					if (Type == LadderType.Parallel)
						return 3 * End1 + (LastHBondDirection == HBondDirection.From0To1 ? 2 : 0);
					else
						return 3 * End1 + (FirstHBondDirection == HBondDirection.From0To1 ? 2 : 0); 
				}
			}

			private BetaLadder(LadderType type, int start0, int end0, int start1, int end1, HBondDirection firstHBondDirection, HBondDirection lastHBondDirection){
				Type=type;
				Start0 = start0;
				End0 = end0;
				Start1 = start1;
				End1 = end1;
				FirstHBondDirection = firstHBondDirection;
				LastHBondDirection = lastHBondDirection;
			}

			public BetaLadder(LadderType type, int res0, int res1, HBondDirection hBondDirection) 
				: this(type, res0, res0, res1, res1, hBondDirection, hBondDirection) {}

			public void AddOneHBond(){
				if (Type == LadderType.Antiparallel) {
					if (LastHBondDirection == HBondDirection.From0To1) {
						LastHBondDirection = HBondDirection.From1To0;
					} else {
						End0 += 2;
						Start1 -= 2;
						LastHBondDirection = HBondDirection.From0To1;
					}
				} else {
					if (LastHBondDirection == HBondDirection.From0To1) {
						End1 += 2;
						LastHBondDirection = HBondDirection.From1To0;
					} else {
						End0 += 2;
						LastHBondDirection = HBondDirection.From0To1;
					}
				}
			}

			public BetaLadder Inverted(){
				HBondDirection invFirstDir = (Type == LadderType.Antiparallel) ? InvertedDirection (LastHBondDirection) : InvertedDirection (FirstHBondDirection);
				HBondDirection invLastDir = (Type == LadderType.Antiparallel) ? InvertedDirection (FirstHBondDirection) : InvertedDirection (LastHBondDirection);
				return new BetaLadder (Type, Start1, End1, Start0, End0, invFirstDir,invLastDir);
			}

		}

		private class BetaBulge {
			public enum BulgeType { Classic, Wide, Antiparallel22, Antiparallel33, Antiparallel15, Antiparallel23, Parallel14, Parallel32, Parallel13, Parallel33 } // New types might be added
			public BulgeType Type {get;private set;}
			public int StartShort{ get; private set; }
			public int EndShort{ get; private set; }
			public int StartLong{ get; private set; }
			public int EndLong{ get; private set; }
			public BetaBulge (BulgeType type, int startShort, int endShort, int startLong, int endLong){
				Type=type;
				StartShort=startShort;
				EndShort=endShort;
				StartLong=startLong;
				EndLong=endLong;
			}
		}

		private class BetaStrandInSheet{
			public SSE SSE { get; set; }
			public int SheetId { get; set; }
			//public int Level { get; set; }
			public bool EvenUp { get; set; }
			public List<BetaStrandInSheet> UpNeighbours { get; set; }
			public List<BetaStrandInSheet> DownNeighbours { get; set; }
			// The ladders must be stored so that strand 0 belongs to this!
			//TODO a check of this condition should be included in AddUpLadder/AddDownLadder
			public List<BetaLadder> UpLadders { get; set; } 
			public List<BetaLadder> DownLadders{ get; set; }
			private bool DFSFlag { get; set; }

			public BetaStrandInSheet(SSE sse, int sheetId, /*int level, */bool evenUp){
				SSE=sse;
				SheetId=sheetId;
				//Level=level;
				EvenUp=evenUp;
				UpNeighbours=new List<BetaStrandInSheet>();
				DownNeighbours=new List<BetaStrandInSheet>();
				UpLadders=new List<BetaLadder>();
				DownLadders=new List<BetaLadder>();
				DFSFlag=false;
			}

			private void DFSAux(bool flagDiscovered, Action<BetaStrandInSheet> actionOnDiscover, Action<BetaStrandInSheet> actionOnReturn){
				if (DFSFlag != flagDiscovered) {
					DFSFlag = flagDiscovered;
					//if (flagDiscovered) Lib.WriteDebug ("[{0}], ", this.SSE.Label);
					actionOnDiscover (this);
					foreach (BetaStrandInSheet s in UpNeighbours.Union (DownNeighbours))
						s.DFSAux (flagDiscovered, actionOnDiscover, actionOnReturn);
					actionOnReturn (this);
				}
			}

			public void DFS (Action<BetaStrandInSheet> actionOnDiscover, Action<BetaStrandInSheet> actionOnReturn){
						//Lib.WriteLineDebug ("DFS from {0}:", this.SSE.Label);
				DFSAux (true, actionOnDiscover, actionOnReturn);
				DFSAux (false, s => {}, s => {});
				//Lib.WriteLineDebug ("");
			}

			public void DFS (Action<BetaStrandInSheet> actionOnDiscover){
				DFS (actionOnDiscover, s => {});
			}

			public override String ToString(){
				return SSE.Label + "[" + SSE.ChainID + SSE.Start + "-" + SSE.End + "]";
			}
		}

		public class HBondSecStrAssigner : ISecStrAssigner{
			public const int MIN_HBONDS_PER_LADDER = 2;
			public const int MIN_HBONDS_PER_HELIX = 1;
			public const int MIN_OVERLAP_FOR_JOINING = 0; // condition for joining 2 beta-strands: end0 >= start1 + MIN_OVERLAP et vice versa
			public const bool STRANDS_BY_ALPHA = true; // if true then residues are assigned to a strand if they have C-alpha included in a cycle (DSSP-style), if false then residues are assigned to a strand if they have any atom included in a cycle (HERA-style)
			public const bool BULGES_BY_ALPHA = false;
			public const bool HELICES_BY_ALPHA = true;
			// To use DSSP-style joining, set MIN_Z_OVERLAP_FOR_JOINING = 3*MIN_OVERLAP_FOR_JOINING+2
			public const int MIN_Z_OVERLAP_FOR_JOINING = 1; // 3*MIN_OVERLAP_FOR_JOINING+2; // condition for joining 2 beta-strands: Z(end0) >= Z(start1) + MIN_OVERLAP et vice versa
			public const bool ALLOW_BULGE_A33 = true;  // antiparallel beta-bulge defined as only 1 missing H-bond from regular beta-ladder (in 2qad chain B ~ resi 15)

			public bool DetectSheets { get; set; }
			public bool DetectHelices { get; set; }

			private List<Residue> residues;
			List<int>[] ignoreHBondsTo_Antiparallel;
			List<int>[] ignoreHBondsFrom_Antiparallel;
			List<int>[] ignoreHBondsTo_Parallel;
			List<int>[] ignoreHBondsFrom_Parallel;

			private HydrogenAdders.IHydrogenAdder hydrogenAdder;
			private IHBondFinder hBondFinder;

			public HBondSecStrAssigner(IEnumerable<Residue> residues, double hBondEnergyCutoff){
				this.DetectSheets=true;
				this.DetectHelices=true;
				this.hydrogenAdder = new HydrogenAdders.DsspLikeAmideHydrogenAdder();
				//List<Atom> atoms = residues.Where (r => r.IsProteinResidueWithCAlpha()).SelectMany (r => r.GetAtoms ()).ToList ();
				//Protein p = hydrogenAdder.AddHydrogens (atoms);
				//this.residues = p.GetResidues ().OrderBy (r=>r).ToList ();
				var validResidues = residues.Where (r => r.HasCAlpha());
				this.residues = hydrogenAdder.AddHydrogens(validResidues).OrderBy (r => r).ToList ();
				ignoreHBondsTo_Antiparallel=new List<int>[residues.Count ()];
				ignoreHBondsFrom_Antiparallel=new List<int>[residues.Count ()];
				ignoreHBondsTo_Parallel=new List<int>[residues.Count ()];
				ignoreHBondsFrom_Parallel=new List<int>[residues.Count ()];
				//this.hBondFinder = new SimpleHBondFinder(this.residues, hBondEnergyCutoff);
				this.hBondFinder = new BoxingHBondFinder(this.residues, hBondEnergyCutoff);
			}

			public HBondSecStrAssigner(Protein p, double hBondEnergyCutoff) : this(p.GetResidues (),hBondEnergyCutoff) {}

			public HBondSecStrAssigner(Chain c, double hBondEnergyCutoff) : this(c.GetResidues (),hBondEnergyCutoff) {}


			private String Ladder2String(BetaLadder ladder){
				char fd = ladder.FirstHBondDirection == BetaLadder.HBondDirection.From0To1 ? 'v' : '^';
				char ld = ladder.LastHBondDirection == BetaLadder.HBondDirection.From0To1 ? 'v' : '^';
				return String.Format ("[{0} {1}-{2} : {3} {4}-{5} {6}{7}]", residues [ladder.Start0].ChainId, residues [ladder.Start0].SeqNumber, residues [ladder.End0].SeqNumber, residues [ladder.Start1].ChainId, residues [ladder.Start1].SeqNumber, residues [ladder.End1].SeqNumber,fd,ld);
			}

			private int Residue2After(int i){
				return ResidueXAfter (i, 2);
			}

			private int Residue2Before(int i){
				return ResidueXBefore (i, 2);
			}

			private int ResidueXAfter(int resIndex, int x){
				if (x < 0)
					return ResidueXBefore (resIndex, -x);
				for (int i = x; i >= 0; i--) {
					if (resIndex + i < residues.Count && residues [resIndex + i].ChainId == residues [resIndex].ChainId && residues [resIndex + i].SeqNumber == residues [resIndex].SeqNumber + x)
						return resIndex + i;
				}
				return -1; //does not exist
			}

			private int ResidueXBefore(int resIndex, int x){
				if (x < 0)
					return ResidueXAfter (resIndex, -x);
				for (int i = x; i >= 0; i--) {
					if (resIndex - i >=0 && residues [resIndex - i].ChainId == residues [resIndex].ChainId && residues [resIndex - i].SeqNumber == residues [resIndex].SeqNumber - x)
						return resIndex - i;
				}
				return -1; //does not exist
			}

			private bool CanExtendLadder(BetaLadder ladder){
				if (ladder.Type == BetaLadder.LadderType.Antiparallel) {
					// Antiparallel
					if (ladder.LastHBondDirection == BetaLadder.HBondDirection.From0To1) {
						return hBondFinder.IsHBond (ladder.Start1, ladder.End0) && (ladder.Start1 > ResidueXAfter (ladder.End0,2));
					} else {
						int i = Residue2After (ladder.End0);
						int j = Residue2Before (ladder.Start1);
						return i>=0 && j>=0 && hBondFinder.IsHBond (i,j);
					}
				} else {
					// Parallel
					if (ladder.LastHBondDirection == BetaLadder.HBondDirection.From0To1) {
						int j = Residue2After (ladder.End1);
						return j >= 0 && hBondFinder.IsHBond (j, ladder.End0);
					} else {
						int i = Residue2After (ladder.End0);
						return i >= 0 && hBondFinder.IsHBond (i, ladder.End1);
					}
				}
			}

			private void AddSafely<T>(ref List<T> list, T newElement){
				if (list == null)
					list = new List<T> ();
				list.Add (newElement);
			}

			private void IgnoreRedundantMotifs(BetaLadder ladder){
				if (ladder.Type == BetaLadder.LadderType.Antiparallel) {
					if (ladder.FirstHBondDirection == BetaLadder.HBondDirection.From0To1)
						AddSafely (ref ignoreHBondsTo_Antiparallel [ladder.Start0], ladder.End1);
					for (int i = 2; i <= ladder.End0 - ladder.Start0; i += 2)
						AddSafely (ref ignoreHBondsTo_Antiparallel [ladder.Start0 + i], ladder.End1 - i);
					for (int i = 0; i < ladder.End0 - ladder.Start0; i += 2)
						AddSafely (ref ignoreHBondsFrom_Antiparallel [ladder.Start0 + i], ladder.End1 - i);
					if (ladder.LastHBondDirection == BetaLadder.HBondDirection.From1To0)
						AddSafely (ref ignoreHBondsFrom_Antiparallel [ladder.End0], ladder.Start1);
				} else { //Parallel
					if (ladder.FirstHBondDirection == BetaLadder.HBondDirection.From0To1) {
						for (int i = 0; i <=ladder.End0-ladder.Start0 ; i += 2) {
							AddSafely (ref ignoreHBondsTo_Parallel [ladder.Start0 + i], ladder.Start1 + i);
						}
						for (int i = 0; i <ladder.End0-ladder.Start0 ; i += 2) {
							AddSafely (ref ignoreHBondsFrom_Parallel [ladder.Start0 + i], ladder.Start1 + 2 + i);
						}
						if (ladder.LastHBondDirection == BetaLadder.HBondDirection.From1To0)
							AddSafely (ref ignoreHBondsFrom_Parallel [ladder.End0], ladder.End1);
					} else {
						for (int i = 2; i <=ladder.End0-ladder.Start0 ; i += 2) {
							AddSafely (ref ignoreHBondsTo_Parallel [ladder.Start0 + i], ladder.Start1 - 2 + i);
						}
						for (int i = 0; i <ladder.End0-ladder.Start0 ; i += 2) {
							AddSafely (ref ignoreHBondsFrom_Parallel [ladder.Start0 + i], ladder.Start1 + i);
						}
						if (ladder.LastHBondDirection == BetaLadder.HBondDirection.From1To0)
							AddSafely (ref ignoreHBondsFrom_Parallel [ladder.End0], ladder.End1);
					}
				}
			}

			private int CheckLadderAndCountHBonds(BetaLadder ladder){
				if (ladder.Type == BetaLadder.LadderType.Antiparallel) {
					if (ladder.End0 - ladder.Start0 != ladder.End1 - ladder.Start1)
						throw new SecStrAssignmentException ("Antiparallel beta-ladder: Strands have different lengths.");
					if ((ladder.End0 - ladder.Start0) % 2 != 0)
						throw new SecStrAssignmentException ("Antiparallel beta-ladder: Strands have even lengths.");
					int numHBonds = (ladder.End0 - ladder.Start0 + 2)
					                - (ladder.FirstHBondDirection == BetaLadder.HBondDirection.From1To0 ? 1 : 0)
					                - (ladder.LastHBondDirection == BetaLadder.HBondDirection.From0To1 ? 1 : 0);
					if (numHBonds<2)
						throw new SecStrAssignmentException ("Antiparallel beta-ladder: Less than 2 hydrogen bonds.");
					return numHBonds;
				} else { //Parallel
					if (ladder.End0 - ladder.Start0 - ladder.End1 + ladder.Start1
					    != (ladder.FirstHBondDirection == BetaLadder.HBondDirection.From1To0 ? 2 : 0) - (ladder.LastHBondDirection == BetaLadder.HBondDirection.From1To0 ? 2 : 0))
						//Lib.WriteLineDebug ("parallel[{0}-{1}|{2}-{3}|{4},{5}]",residues[ladder.Start0],residues[ladder.End0],residues[ladder.Start1],residues[ladder.End1],ladder.FirstHBondDirection,ladder.LastHBondDirection);
						throw new SecStrAssignmentException ("Parallel beta-ladder: Strands have incompatible lengths.");
					int numHBonds = 1 + (ladder.End0 - ladder.Start0 + ladder.End1 - ladder.Start1) / 2;
					if (numHBonds<2)
						throw new SecStrAssignmentException ("Parallel beta-ladder: Less than 2 hydrogen bonds.");
					return numHBonds;
				}
			}

			private bool DoOverlap(BetaStrandInSheet a, BetaStrandInSheet b){
				if (a.SSE.ChainID == b.SSE.ChainID) {
					foreach (BetaLadder la in a.UpLadders.Union (a.DownLadders))
						foreach (BetaLadder lb in b.UpLadders.Union (b.DownLadders))
							if (la.ZEnd0 - lb.ZStart0 >= MIN_Z_OVERLAP_FOR_JOINING && lb.ZEnd0 - la.ZStart0 >= MIN_Z_OVERLAP_FOR_JOINING)
								return true;
				}
				return false;
				//return a.SSE.ChainID==b.SSE.ChainID && a.SSE.End >= b.SSE.Start + MIN_OVERLAP_FOR_JOINING && b.SSE.End >= a.SSE.Start + MIN_OVERLAP_FOR_JOINING;
			}

			private bool AreShortSideOfBetaBulge(BetaStrandInSheet a, BetaStrandInSheet b){
				if (a.SSE.ChainID != b.SSE.ChainID)
					return false;
				foreach (BetaLadder la in a.UpLadders.Union (a.DownLadders)) {
					foreach (BetaLadder lb in b.UpLadders.Union (b.DownLadders)) {
						if (BuildBetaBulge (la, lb) != null)
							return true;
					}
				}
				return false;
			}

			private bool AreLongSideOfBetaBulge(BetaStrandInSheet a, BetaStrandInSheet b){
				if (a.SSE.ChainID != b.SSE.ChainID)
					return false;
				foreach (BetaLadder la in a.UpLadders.Union (a.DownLadders)) {
					foreach (BetaLadder lb in b.UpLadders.Union (b.DownLadders)) {
						if (BuildBetaBulge (la.Inverted (), lb.Inverted ()) != null)
							return true;
					}
				}
				return false;
			}

			/* Tries to build a beta-bulge so than the la.Strand0 and lb.Strand0 form the SHORT side of the bulge. */
			private BetaBulge BuildBetaBulge(BetaLadder la,BetaLadder lb){
				return BuildBetaBulgeWithFirstLadderFirst (la, lb) ?? BuildBetaBulgeWithFirstLadderFirst (lb, la);
			}
			/* Tries to build a beta-bulge so than the la.Strand0 and lb.Strand0 form the SHORT side of the bulge. */
			private BetaBulge BuildBetaBulgeWithFirstLadderFirst(BetaLadder la,BetaLadder lb){
				if (residues [la.Start1].ChainId == residues [lb.Start1].ChainId) {
					if (la.Type == BetaLadder.LadderType.Antiparallel && lb.Type == BetaLadder.LadderType.Antiparallel) {
						//Antiparallel types
						if (lb.Start0 == la.End0 && la.Start1 == ResidueXAfter (lb.End1, 1)) {
							// classic bulge a-b (type A12)
							return new BetaBulge (BetaBulge.BulgeType.Classic, la.End0, lb.Start0, lb.End1, la.Start1);
						}
						if (lb.Start0 == ResidueXAfter (la.End0, 2) && la.Start1 == ResidueXAfter (lb.End1, 3)
							&& la.LastHBondDirection == BetaLadder.HBondDirection.From1To0 && lb.FirstHBondDirection == BetaLadder.HBondDirection.From0To1) {
							// wide bulge a-b (type A34)
							return new BetaBulge (BetaBulge.BulgeType.Wide, la.End0, lb.Start0, lb.End1, la.Start1);
						}
						if (lb.Start0 == ResidueXAfter (la.End0, 1) && la.Start1 == ResidueXAfter (lb.End1, 1)
							&& la.LastHBondDirection == BetaLadder.HBondDirection.From0To1 && lb.FirstHBondDirection == BetaLadder.HBondDirection.From0To1) {
							// bulge-not-bulge a-b (type A22) found in 1gei
							return new BetaBulge (BetaBulge.BulgeType.Antiparallel22, la.End0, lb.Start0, lb.End1, la.Start1);
						}
						if (ALLOW_BULGE_A33 && lb.Start0 == ResidueXAfter (la.End0, 2) && la.Start1 == ResidueXAfter (lb.End1, 2)
							&& la.LastHBondDirection == BetaLadder.HBondDirection.From0To1 && lb.FirstHBondDirection == BetaLadder.HBondDirection.From0To1) {
							// bulge-not-bulge a-b (type A33)
							return new BetaBulge (BetaBulge.BulgeType.Antiparallel33, la.End0, lb.Start0, lb.End1, la.Start1);
						}
						if (lb.Start0 == la.End0 && la.Start1 == ResidueXAfter (lb.End1, 4)
							&& la.LastHBondDirection == BetaLadder.HBondDirection.From1To0 && lb.FirstHBondDirection == BetaLadder.HBondDirection.From1To0) {
							// arch-bulge a-b (type A15) found in 1gjm
							return new BetaBulge (BetaBulge.BulgeType.Antiparallel15, la.End0, lb.Start0, lb.End1, la.Start1);
						}
						if (lb.Start0 == ResidueXAfter (la.End0, 1) && la.Start1 == ResidueXAfter (lb.End1, 2)
							&& la.LastHBondDirection == BetaLadder.HBondDirection.From0To1 && lb.FirstHBondDirection == BetaLadder.HBondDirection.From0To1) {
							// arch-bulge a-b (type A15) found in 1gjm
							return new BetaBulge (BetaBulge.BulgeType.Antiparallel23, la.End0, lb.Start0, lb.End1, la.Start1);
						}
					}
					if (la.Type == BetaLadder.LadderType.Parallel && lb.Type == BetaLadder.LadderType.Parallel) {
						//Parallel types
						if (lb.Start0 == la.End0 && lb.Start1 == ResidueXAfter (la.End1, 3)
							&& la.LastHBondDirection == BetaLadder.HBondDirection.From0To1 && lb.FirstHBondDirection == BetaLadder.HBondDirection.From1To0) {
							// parallel bulge (type P14)
							return new BetaBulge (BetaBulge.BulgeType.Parallel14, la.End0, lb.Start0, la.End1, lb.Start1);
						}
						if (lb.Start0 == ResidueXAfter (la.End0, 2) && lb.Start1 == ResidueXAfter (la.End1, 1)
							&& la.LastHBondDirection == BetaLadder.HBondDirection.From1To0 && lb.FirstHBondDirection == BetaLadder.HBondDirection.From0To1) {
							// parallel bulge (type P32)
							return new BetaBulge (BetaBulge.BulgeType.Parallel32, la.End0, lb.Start0, la.End1, lb.Start1);
						}
						if (lb.Start0 == la.End0 && lb.Start1 == ResidueXAfter (la.End1, 2)
							&& la.LastHBondDirection == BetaLadder.HBondDirection.From0To1 && lb.FirstHBondDirection == BetaLadder.HBondDirection.From0To1) {
							// parallel bulge (type P13)
							return new BetaBulge (BetaBulge.BulgeType.Parallel13, la.End0, lb.Start0, la.End1, lb.Start1);
						}
						if (lb.Start0 == ResidueXAfter (la.End0, 2) && lb.Start1 == ResidueXAfter (la.End1, 2)
							&& la.LastHBondDirection == BetaLadder.HBondDirection.From1To0 && lb.FirstHBondDirection == BetaLadder.HBondDirection.From0To1) {
							// parallel bulge (type P33)
							return new BetaBulge (BetaBulge.BulgeType.Parallel33, la.End0, lb.Start0, la.End1, lb.Start1);
						}
					}
				}
				return null;
			}

			private void InvertSheetLevels(BetaStrandInSheet strand){
				Lib.WriteLineDebug ("Invert sheet: {0}.",strand.SheetId);
				strand.DFS (s => {
					//s.Level = -s.Level;
					s.EvenUp = !s.EvenUp;
					var aux = s.DownNeighbours;
					s.DownNeighbours = s.UpNeighbours;
					s.UpNeighbours = aux;
					var aux2 = s.DownLadders;
					s.DownLadders = s.UpLadders;
					s.UpLadders = aux2;
				});
			}

			private void SetIdToSheet(BetaStrandInSheet strand, int newID){
				Lib.WriteLineDebug ("Set ID: {0} to {1}.",strand.SheetId,newID);
				strand.DFS (s => {
					s.SheetId = newID;
				});
				Lib.WriteLineDebug ("");
			}

			private void LinkStrands(BetaStrandInSheet lower, BetaStrandInSheet upper, BetaLadder ladder){
				Lib.WriteLineDebug ("Linking strands {0} ({1}) and {2} ({3}) {4}.",lower.SSE.Label, lower.SheetId,upper.SSE.Label, upper.SheetId, Ladder2String (ladder));
				if (lower == upper)
					Lib.WriteWarning ("Trying to link a beta-strand to itself.");
				if (lower.SheetId != upper.SheetId)
					throw new ArgumentException (); 
				lower.UpNeighbours.Add (upper);
				lower.UpLadders.Add (ladder);
				upper.DownNeighbours.Add (lower);
				upper.DownLadders.Add (ladder.Inverted());
			}

			private void UnlinkStrands(BetaStrandInSheet lower, BetaStrandInSheet upper){
				if (lower == upper)
					Lib.WriteWarning ("Trying to unlink a beta-strand from itself.");
				int i = lower.UpNeighbours.IndexOf (upper);
				int j = upper.DownNeighbours.IndexOf (lower);
				Lib.WriteLineDebug ("Unlinking: upper={0}, lower={1}", i, j);
				lower.UpNeighbours.RemoveAt (i);
				lower.UpLadders.RemoveAt (i);
				upper.DownNeighbours.RemoveAt (j);
				upper.DownLadders.RemoveAt (j);
			}

			private void AddVertexAndPossiblyMerge(List<BetaStrandInSheet> vertices, BetaStrandInSheet newVertex){
				Lib.WriteLineDebug ("Adding vertex {0}.",newVertex);
				vertices.Add (newVertex);
				int iNew = vertices.Count - 1;
				for (int i = vertices.Count - 2; i >= 0; i--) {
					if (DoOverlap (vertices [i], vertices [iNew]) || AreShortSideOfBetaBulge (vertices [i], vertices [iNew]) || AreLongSideOfBetaBulge (vertices [i], vertices [iNew])) {
						bool merged = TryMergeSheetsByStrands (vertices [i], vertices [iNew]);
						if (merged) {
							vertices.RemoveAt (iNew);
							iNew = i;
						}		
					}
				}
			}

			private bool TryMergeSheetsByStrands(BetaStrandInSheet extendedStrand, BetaStrandInSheet extendingStrand){
				Lib.WriteLineDebug ("Merging strands {0} ({1}) and {2} ({3}).", extendedStrand.SSE.Label, extendedStrand.SheetId, extendingStrand.SSE.Label, extendingStrand.SheetId);
				if (extendedStrand == extendingStrand)
					throw new ArgumentException ();
				string chain = extendedStrand.SSE.ChainID;
				int start = Math.Min (extendedStrand.SSE.Start, extendingStrand.SSE.Start);
				int end = Math.Max (extendedStrand.SSE.End, extendingStrand.SSE.End);

				if (extendedStrand.SheetId == extendingStrand.SheetId) {
					// in the same sheet
					if (extendedStrand.EvenUp != extendingStrand.EvenUp) {
						Lib.WriteLineDebug ("Important! Cycle with inconsistent direction would arise in beta-strand graph (around chain {0} {1}-{2})! Skipping merging strands.", chain, start, end);
						//Lib.WriteWarning ("Important! Cycle with inconsistent direction would arise in beta-strand graph (around chain {0} {1}-{2})! Skipping merging strands.", chain, start, end);
						//return false;
					} else {
						Lib.WriteLineDebug ("Cycle detected in beta-strand graph (around chain {0} {1}-{2}).", chain, start, end);
					}
					/*if (extendedStrand.Level != extendingStrand.Level) {
						Lib.WriteWarning ("Cycle with inconsistent level (might be a barrel) detected in beta-strand graph (around chain {0} {1}-{2})! Skipping merging strands.", extendedStrand.SSE.ChainID, Math.Min (extendedStrand.SSE.Start, extendingStrand.SSE.Start), Math.Max (extendedStrand.SSE.End, extendingStrand.SSE.End));
					}*/
				} else {
					// in different sheets
					if (extendedStrand.EvenUp != extendingStrand.EvenUp)
						InvertSheetLevels (extendingStrand);
					SetIdToSheet (extendingStrand, extendedStrand.SheetId);
				}

				extendedStrand.SSE.Start = start;
				extendedStrand.SSE.End = end;
				BetaStrandInSheet[] oldUpNeighbours = extendingStrand.UpNeighbours.ToArray ();
				BetaLadder[] oldUpLadders = extendingStrand.UpLadders.ToArray ();
				BetaStrandInSheet[] oldDownNeighbours = extendingStrand.DownNeighbours.ToArray ();
				BetaLadder[] oldDownLadders = extendingStrand.DownLadders.ToArray ();
				Lib.WriteLineDebug ("Extended: {0}",extendedStrand.SSE);
				Lib.WriteLineDebug ("Extending: {0}",extendingStrand.SSE);
				for (int i = 0; i < oldUpNeighbours.Length; i++) {
					UnlinkStrands (extendingStrand, oldUpNeighbours [i]);
					LinkStrands (extendedStrand, oldUpNeighbours [i], oldUpLadders [i]);
				}
				for (int i = 0; i < oldDownNeighbours.Length; i++) {
					UnlinkStrands (oldDownNeighbours [i], extendingStrand);
					LinkStrands (oldDownNeighbours [i], extendedStrand, oldDownLadders [i].Inverted ());
				}
				return true;
			}

			private SSE GetStrand0(BetaLadder ladder){
				string chainId = residues [ladder.Start0].ChainId;
				int start = residues [ladder.Start0].SeqNumber;
				int end = residues [ladder.End0].SeqNumber;
				if (STRANDS_BY_ALPHA) {
					if (ladder.FirstHBondDirection == BetaLadder.HBondDirection.From1To0)
						start++;
					if (ladder.LastHBondDirection == BetaLadder.HBondDirection.From0To1)
						end--;
				}
				return new SSE (null, chainId, start, end, LadderSSEType (ladder), null);
			}

			private SSE GetStrand1(BetaLadder ladder){
				string chainId = residues [ladder.Start1].ChainId;
				int start = residues [ladder.Start1].SeqNumber;
				int end = residues [ladder.End1].SeqNumber;
				if (STRANDS_BY_ALPHA) {
					if (ladder.Type == BetaLadder.LadderType.Antiparallel) {
						if (ladder.FirstHBondDirection == BetaLadder.HBondDirection.From1To0)
							end--;
						if (ladder.LastHBondDirection == BetaLadder.HBondDirection.From0To1)
							start++;
					}
					if (ladder.Type == BetaLadder.LadderType.Parallel) {
						if (ladder.FirstHBondDirection == BetaLadder.HBondDirection.From0To1)
							start++;
						if (ladder.LastHBondDirection == BetaLadder.HBondDirection.From1To0)
							end--;
					}
				}
				return new SSE (null, chainId, start, end, LadderSSEType (ladder), null);
			}

			private SSE GetShortStrand(BetaBulge bulge){
				if (BULGES_BY_ALPHA) {
					throw new NotImplementedException ();
				} else {
					switch (bulge.Type) {
					case BetaBulge.BulgeType.Classic:
						return new SSE (null, residues [bulge.StartShort].ChainId, residues [bulge.StartShort].SeqNumber, residues [bulge.EndShort].SeqNumber, SSE.BULGE_CLASSIC_SHORT_SIDE_TYPE, null);
					case BetaBulge.BulgeType.Wide:
						return new SSE (null, residues [bulge.StartShort].ChainId, residues [bulge.StartShort].SeqNumber, residues [bulge.EndShort].SeqNumber, SSE.BULGE_WIDE_SHORT_SIDE_TYPE, null);
					case BetaBulge.BulgeType.Antiparallel22:
						return new SSE (null, residues [bulge.StartShort].ChainId, residues [bulge.StartShort].SeqNumber, residues [bulge.EndShort].SeqNumber, SSE.BULGE_ANTIPARALLEL22_SHORT_SIDE_TYPE, null);
					case BetaBulge.BulgeType.Antiparallel33:
						return new SSE (null, residues [bulge.StartShort].ChainId, residues [bulge.StartShort].SeqNumber, residues [bulge.EndShort].SeqNumber, SSE.BULGE_ANTIPARALLEL33_SHORT_SIDE_TYPE, null);
					case BetaBulge.BulgeType.Antiparallel15:
						return new SSE (null, residues [bulge.StartShort].ChainId, residues [bulge.StartShort].SeqNumber, residues [bulge.EndShort].SeqNumber, SSE.BULGE_ANTIPARALLEL15_SHORT_SIDE_TYPE, null);
					case BetaBulge.BulgeType.Antiparallel23:
						return new SSE (null, residues [bulge.StartShort].ChainId, residues [bulge.StartShort].SeqNumber, residues [bulge.EndShort].SeqNumber, SSE.BULGE_ANTIPARALLEL23_SHORT_SIDE_TYPE, null);
					case BetaBulge.BulgeType.Parallel14:
						return new SSE (null, residues [bulge.StartShort].ChainId, residues [bulge.StartShort].SeqNumber, residues [bulge.EndShort].SeqNumber, SSE.BULGE_PARALLEL14_SHORT_SIDE_TYPE, null);
					case BetaBulge.BulgeType.Parallel32:
						return new SSE (null, residues [bulge.StartShort].ChainId, residues [bulge.StartShort].SeqNumber, residues [bulge.EndShort].SeqNumber, SSE.BULGE_PARALLEL32_SHORT_SIDE_TYPE, null);
					case BetaBulge.BulgeType.Parallel13:
						return new SSE (null, residues [bulge.StartShort].ChainId, residues [bulge.StartShort].SeqNumber, residues [bulge.EndShort].SeqNumber, SSE.BULGE_PARALLEL13_SHORT_SIDE_TYPE, null);
					case BetaBulge.BulgeType.Parallel33:
						return new SSE (null, residues [bulge.StartShort].ChainId, residues [bulge.StartShort].SeqNumber, residues [bulge.EndShort].SeqNumber, SSE.BULGE_PARALLEL33_SHORT_SIDE_TYPE, null);
					default:
						throw new NotImplementedException ();
					}
				}
			}
		
			private SSE GetLongStrand(BetaBulge bulge){
				if (BULGES_BY_ALPHA) {
					throw new NotImplementedException ();
				} else {
					switch (bulge.Type) {
					case BetaBulge.BulgeType.Classic:
						return new SSE (null, residues [bulge.StartLong].ChainId, residues [bulge.StartLong].SeqNumber, residues [bulge.EndLong].SeqNumber, SSE.BULGE_CLASSIC_LONG_SIDE_TYPE, null);
					case BetaBulge.BulgeType.Wide:
						return new SSE (null, residues [bulge.StartLong].ChainId, residues [bulge.StartLong].SeqNumber, residues [bulge.EndLong].SeqNumber, SSE.BULGE_WIDE_LONG_SIDE_TYPE, null);
					case BetaBulge.BulgeType.Antiparallel22:
						return new SSE (null, residues [bulge.StartLong].ChainId, residues [bulge.StartLong].SeqNumber, residues [bulge.EndLong].SeqNumber, SSE.BULGE_ANTIPARALLEL22_LONG_SIDE_TYPE, null);
					case BetaBulge.BulgeType.Antiparallel33:
						return new SSE (null, residues [bulge.StartLong].ChainId, residues [bulge.StartLong].SeqNumber, residues [bulge.EndLong].SeqNumber, SSE.BULGE_ANTIPARALLEL33_LONG_SIDE_TYPE, null);
					case BetaBulge.BulgeType.Antiparallel15:
						return new SSE (null, residues [bulge.StartLong].ChainId, residues [bulge.StartLong].SeqNumber, residues [bulge.EndLong].SeqNumber, SSE.BULGE_ANTIPARALLEL15_LONG_SIDE_TYPE, null);
					case BetaBulge.BulgeType.Antiparallel23:
						return new SSE (null, residues [bulge.StartLong].ChainId, residues [bulge.StartLong].SeqNumber, residues [bulge.EndLong].SeqNumber, SSE.BULGE_ANTIPARALLEL23_LONG_SIDE_TYPE, null);
					case BetaBulge.BulgeType.Parallel14:
						return new SSE (null, residues [bulge.StartLong].ChainId, residues [bulge.StartLong].SeqNumber, residues [bulge.EndLong].SeqNumber, SSE.BULGE_PARALLEL14_LONG_SIDE_TYPE, null);
					case BetaBulge.BulgeType.Parallel32:
						return new SSE (null, residues [bulge.StartLong].ChainId, residues [bulge.StartLong].SeqNumber, residues [bulge.EndLong].SeqNumber, SSE.BULGE_PARALLEL32_LONG_SIDE_TYPE, null);
					case BetaBulge.BulgeType.Parallel13:
						return new SSE (null, residues [bulge.StartLong].ChainId, residues [bulge.StartLong].SeqNumber, residues [bulge.EndLong].SeqNumber, SSE.BULGE_PARALLEL13_LONG_SIDE_TYPE, null);
					case BetaBulge.BulgeType.Parallel33:
						return new SSE (null, residues [bulge.StartLong].ChainId, residues [bulge.StartLong].SeqNumber, residues [bulge.EndLong].SeqNumber, SSE.BULGE_PARALLEL33_LONG_SIDE_TYPE, null);
					default:
						throw new NotImplementedException ();
					}
				}
			}

			private char LadderSSEType(BetaLadder ladder){
				if (ladder.Start0 == ladder.Start1) {
					// C7 motif
					if (residues [ladder.End0].SeqNumber - residues [ladder.Start0].SeqNumber != 2)
						throw new  SecStrAssignmentException ("C7 turn with more than 1 stabilizing H-bond.");
					else
						return SSE.TURN_C7_TYPE;
				} else {
					// real beta ladder
					return (CheckLadderAndCountHBonds (ladder) > 2) ? SSE.SHEET_TYPE : SSE.ISOLATED_BETA_BRIDGE_TYPE;
				}
			}

			public SecStrAssignment GetSecStrAssignment(){
				if (DetectSheets) {
					if (DetectHelices) {
						return SecStrAssignment.Combine (GetSheets (), GetHelices ());
					} else {
						return GetSheets ();
					}
				} else {
					if (DetectHelices) {
						return GetHelices ();
					} else {
						return new SecStrAssignment (new SSE[0]);
					}
				}
			}

			private SecStrAssignment GetHelices(){
				var gHelices = GetXHelices (3, SSE.HELIX_G_TYPE);
				var hHelices = GetXHelices (4, SSE.HELIX_H_TYPE);
				var iHelices = GetXHelices (5, SSE.HELIX_I_TYPE);
				return new SecStrAssignment (gHelices.Concat (hHelices).Concat (iHelices));
			}

			/** Returns the list of helices formed by (i+x -> i) H bonds. */
			private List<SSE> GetXHelices(int x, char assignedType){
				List<SSE> helices = new List<SSE> ();
				int currentStart = -1; // -1 = currently not making helix
				int currentEnd = -1;
				for (int i = 0; i < residues.Count; i++) {
					int j = ResidueXAfter (i, x);
					if (currentStart == -1 && hBondFinder.IsHBond (j, i)) { 
						//start helix
						currentStart = i;
						currentEnd = j;
					} else if (currentStart >= 0 && hBondFinder.IsHBond (j, i)) { 
						//continue helix
						currentEnd = j;
					} else if (currentStart >= 0 && !hBondFinder.IsHBond (j, i)) {
						//finish or reject helix
						if (i - currentStart >= MIN_HBONDS_PER_HELIX) {
							helices.Add (new SSE (null, 
								residues [currentStart].ChainId, 
								residues [currentStart].SeqNumber + (HELICES_BY_ALPHA ? 1 : 0), 
								residues [currentEnd].SeqNumber - (HELICES_BY_ALPHA ? 1 : 0), 
								assignedType, null));
						}
						currentStart = -1;
					}
				}
				return helices;
			}

			private SecStrAssignment GetSheets(){
				// Find H-binding motifs
				DateTime t0 = DateTime.Now;
				List<BetaLadder> ladders = new List<BetaLadder> ();
				List<Tuple<Residue,Residue>> hBonds = new List<Tuple<Residue, Residue>> ();
				for (int i = 0; i < residues.Count; i++) {
					int nLaddersBefore = ladders.Count;
					List<int> hAcceptors = hBondFinder.FindHAcceptors (i).Where (j => j > i)/*.Where (j => ignoreHBondsTo [i] == null || !ignoreHBondsTo [i].Contains (j))*/.ToList ();
					foreach (int j in hAcceptors) {
						// try motif A
						if ((ignoreHBondsTo_Antiparallel [i] == null || !ignoreHBondsTo_Antiparallel [i].Contains (j)) && hBondFinder.IsHBond (j, i)) {
							BetaLadder ladder = new BetaLadder (BetaLadder.LadderType.Antiparallel, i, j, BetaLadder.HBondDirection.From0To1);
							ladder.AddOneHBond ();
							while (CanExtendLadder (ladder))
								ladder.AddOneHBond ();
							CheckLadderAndCountHBonds (ladder);
							IgnoreRedundantMotifs (ladder);
							ladders.Add (ladder);
						}						
						// try motif C
						if ((ignoreHBondsTo_Parallel [i] == null || !ignoreHBondsTo_Parallel [i].Contains (j)) && hBondFinder.IsHBond (Residue2After (j), i)) {
							BetaLadder ladder = new BetaLadder (BetaLadder.LadderType.Parallel, i, j, BetaLadder.HBondDirection.From0To1);
							ladder.AddOneHBond ();
							while (CanExtendLadder (ladder))
								ladder.AddOneHBond ();
							CheckLadderAndCountHBonds (ladder);
							IgnoreRedundantMotifs (ladder);
							ladders.Add (ladder);
						}		
					}
					List<int> hDonors = hBondFinder.FindHDonors (i).Where (j => j > i)/*.Where (j => ignoreHBondsFrom [i] == null || !ignoreHBondsFrom [i].Contains (j))*/.ToList ();
					foreach (int j in hDonors) {
						// try motif B
						if ((ignoreHBondsFrom_Antiparallel [i] == null || !ignoreHBondsFrom_Antiparallel [i].Contains (j)) && hBondFinder.IsHBond (Residue2After (i), Residue2Before (j))) {
							BetaLadder ladder = new BetaLadder (BetaLadder.LadderType.Antiparallel, i, j, BetaLadder.HBondDirection.From1To0);
							ladder.AddOneHBond ();
							while (CanExtendLadder (ladder))
								ladder.AddOneHBond ();
							CheckLadderAndCountHBonds (ladder);
							IgnoreRedundantMotifs (ladder);
							ladders.Add (ladder);
						}				
						// try motif D
						if ((ignoreHBondsFrom_Parallel [i] == null || !ignoreHBondsFrom_Parallel [i].Contains (j)) && hBondFinder.IsHBond (Residue2After (i), j)) {
							BetaLadder ladder = new BetaLadder (BetaLadder.LadderType.Parallel, i, j, BetaLadder.HBondDirection.From1To0);
							ladder.AddOneHBond ();
							while (CanExtendLadder (ladder))
								ladder.AddOneHBond ();
							CheckLadderAndCountHBonds (ladder);
							IgnoreRedundantMotifs (ladder);
							ladders.Add (ladder);
						}
					}
					if (ladders.Count - nLaddersBefore > 1)
						Lib.WriteLineDebug ("More than one beta-ladder motif found from residue {0}.", residues [i].ToString (true));
					hBonds.AddRange (hAcceptors.Select (a => new Tuple<Residue,Residue> (residues [i], residues [a])));
					hBonds.AddRange (hDonors.Select (d => new Tuple<Residue,Residue> (residues [d], residues [i])));
				}
				Lib.WriteLineDebug ("Time for finding H-bond patterns: " + DateTime.Now.Subtract (t0));

				foreach (BetaLadder l in ladders)
					if (l.Start0 == l.Start1 && (l.End0 - l.Start0) > 2)
						Lib.WriteWarning ("Strange secondary structure in chain {0} {1}-{2}.", residues [l.Start0].ChainId, residues [l.Start0].SeqNumber, residues [l.End0].SeqNumber);

				List<SSE> c7Turns = ladders
					.Where (l => l.Start0 == l.Start1)
					.OrderBy (l => l.Start0)
					.Select (l => GetStrand0 (l)/*new SSE ("x", residues [l.Start0].ChainID, residues [l.Start0].ResSeq, residues [l.End0].ResSeq, SSE.TURN_C7_TYPE,null)*/)
					.ToList ();

				// Build C7 turns and wiggles
				List<SSE> c7TurnsAndWiggles = new List<SSE> ();
				string currentChain = null;
				int? firstResi = null;
				int? lastResi = null;
				int strandCounter = 0;
				foreach (SSE t in c7Turns) { 
					if (t.ChainID == currentChain && t.End == lastResi + 1) { 
						// elongation of existing SSE
						lastResi++;
					} else {
						if (currentChain != null) { 
							// finishing the old SSE
							char type = (lastResi - firstResi > 2) ? SSE.WIGGLE_C7_TYPE : SSE.TURN_C7_TYPE;
							c7TurnsAndWiggles.Add (new SSE ("" + type + (strandCounter++), currentChain, (int)firstResi, (int)lastResi, type, null));
						}
						// starting a new SSE
						currentChain = t.ChainID;
						firstResi = t.Start;
						lastResi = t.End;
					}
				}

				// DEBUG: Write some info about the ladders
				if (Lib.DoWriteDebug) {
					foreach (BetaLadder ladder in ladders) {
						if (ladder.Start0 != ladder.Start1)
							Lib.WriteInColor (ConsoleColor.Cyan, "*");
						else
							Lib.WriteInColor (ConsoleColor.Red, "*");
						Lib.WriteLineDebug (Ladder2String (ladder));
					}
				}

				// Build beta-sheet graph
				List<BetaLadder> realLadders = ladders.Where (l => l.Start0 != l.Start1).Where (l => CheckLadderAndCountHBonds (l) >= MIN_HBONDS_PER_LADDER).ToList ();
				List<BetaStrandInSheet> vertices = new List<BetaStrandInSheet> ();
				int sheetCounter = 0;
				foreach (BetaLadder ladder in realLadders) {
					SSE strand0 = GetStrand0 (ladder).RelabeledCopy ("x" + strandCounter++);
					SSE strand1 = GetStrand1 (ladder).RelabeledCopy ("x" + strandCounter++);
					BetaStrandInSheet lower = new BetaStrandInSheet (strand0, sheetCounter, /*0, */strand0.Start % 2 == 0);
					BetaStrandInSheet upper = new BetaStrandInSheet (strand1, sheetCounter++, /*1, */strand1.Start % 2 != 0);
					LinkStrands (lower, upper, ladder);
					AddVertexAndPossiblyMerge (vertices, lower);
					AddVertexAndPossiblyMerge (vertices, upper);
				}

				// Get beta-sheets (i.e. connected components of the beta-sheet graph)
				List<BetaStrandInSheet> seeds = new List<BetaStrandInSheet> ();
				foreach (BetaStrandInSheet v in vertices) {
					if (!seeds.Any (u => u.SheetId == v.SheetId)) {
						//SetIdToSheet (v, sheetCounter++);
						seeds.Add (v);
					}
				}
				seeds.Sort ((p, q) => SSE.Compare (p.SSE, q.SSE));
				List<List<BetaStrandInSheet>> sheets = new List<List<BetaStrandInSheet>> ();
				sheetCounter = 1;
				foreach (BetaStrandInSheet v in seeds) {
					List<BetaStrandInSheet> sheet = new List<BetaStrandInSheet> ();
					v.DFS (u => {
						u.SheetId = sheetCounter;
						u.SSE.SheetId = sheetCounter;
						IEnumerable<BetaLadder> uLadders = u.DownLadders.Union (u.UpLadders);
						u.SSE.Type = (u.SSE.End - u.SSE.Start > 0 || uLadders.Any (l => LadderSSEType (l) == SSE.SHEET_TYPE)) ? SSE.SHEET_TYPE : SSE.ISOLATED_BETA_BRIDGE_TYPE;
						sheet.Add (u);
					});
					sheetCounter++;
					sheets.Add (sheet);
				}

				// // DEBUG: Write out sheet info.
				// Lib.WriteLineDebug ("From list of sheets:");
				// foreach (var sheet in sheets) {
				// 	Lib.WriteLineDebug ("Sheet {0}:", sheet.First ().SheetId);
				// 	foreach (var u in sheet) {
				// 		Lib.WriteLineDebug ("Strand {0} ({1} {2}-{3}): up: {4}, down: {5}",
				// 			u.SSE.Label,
				// 			u.SSE.ChainID,
				// 			u.SSE.Start,
				// 			u.SSE.End,
				// 			Lib.EnumerateWithCommas (u.UpNeighbours.Select (n => n.SSE.Label)),
				// 			Lib.EnumerateWithCommas (u.DownNeighbours.Select (n => n.SSE.Label)));
				// 	}
				// }

				List<SSE> realStrands = sheets
					.SelectMany (s => s)
					.Select (v => v.SSE)
					.OrderBy (sse => sse)
					.Select ((sse, i) => sse.RelabeledCopy (sse.Type.ToString () + i))
					.ToList ();
				
				// Get beta-bulges //TODO: Beta-bulges can be detected directly (as H-bond pattern) instead of this way.
				List<BetaBulge> bulges = new List<BetaBulge> ();
				BetaLadder[] bothSideRealLadders = realLadders.Union (realLadders.Select (l => l.Inverted ())).ToArray ();
				for (int i = 0; i < bothSideRealLadders.Length; i++) {
					for (int j = i + 1; j < bothSideRealLadders.Length; j++) {
						BetaBulge bulge = BuildBetaBulge (bothSideRealLadders [i], bothSideRealLadders [j]);
						if (bulge != null)
							bulges.Add (bulge);
					}
				}

				int bulgeCounter = 0;

				//List<SSE> bulgeSides=bulges.SelectMany(b=>new SSE[]{GetShortStrand (b),GetLongStrand (b)}).ToList ();
				List<SSE> bulgeSides = new List<SSE> ();
				//Lib.WriteLineDebug ("Detected bulges:");
				foreach (BetaBulge bulge in bulges) {
					foreach (SSE b in new SSE[]{ GetShortStrand (bulge), GetLongStrand (bulge) }) {
						SSE side = b.RelabeledCopy (b.Type.ToString () + bulgeCounter.ToString ());
						SSE[] includingStrands = realStrands.Where (s => s.ChainID == side.ChainID && s.Start <= side.Start && s.End >= side.End).ToArray ();
						if (includingStrands.Length != 1) {
							throw new InvalidOperationException ("One side of a beta-bulge belongs to more than one beta-strand (this should never happen)!");
						} else {
							side.SheetId = includingStrands [0].SheetId;
							includingStrands [0].AddNestedSSE (side);
							bulgeSides.Add (side);
						}
					}
					bulgeCounter++;
				}

				List<SSE> resultSSEs = c7TurnsAndWiggles.Union (realStrands)/*.Union (bulgeSides)*/.OrderBy (sse => sse).ToList ();
				/*foreach (SSE sse in result) {
					int s = ResidueXBefore (sse.Start, 1);
					int e = ResidueXAfter (sse.End, 1);
					if (s == -1 || e == -1) {
						sse.AddComment ("Edge of chain. ");
					} else if (Enumerable.Range (s, e - s + 1).Any (i => residues [i].GetAtoms ().Any (a => a.AltLoc != ' '))){
						sse.AddComment ("Alternative locations. ");
					}
				}*/

				Dictionary<Tuple<string,int,int>,int> chainStartEnd2Index = new Dictionary<Tuple<string, int, int>, int> ();
				for (int i = 0; i < resultSSEs.Count; i++) {
					SSE sse = resultSSEs [i];
					chainStartEnd2Index [new Tuple<string,int,int> (sse.ChainID, sse.Start, sse.End)] = i;
				}
				List<Tuple<int,int,int>> edges = new List<Tuple<int, int,int>> ();
				foreach (var seed in seeds) {
					seed.DFS (u => {
						int vertex1 = chainStartEnd2Index [new Tuple<string,int,int> (u.SSE.ChainID, u.SSE.Start, u.SSE.End)];
						for (int i = 0; i < u.DownNeighbours.Count; i++) {
							var v = u.DownNeighbours [i];
							var ladder = u.DownLadders [i];
							int vertex2 = chainStartEnd2Index [new Tuple<string,int,int> (v.SSE.ChainID, v.SSE.Start, v.SSE.End)];
							int ladderType = ladder.Type == BetaLadder.LadderType.Parallel ? 1 : -1;
							edges.Add (new Tuple<int,int,int> (Math.Min (vertex1, vertex2), Math.Max (vertex1, vertex2), ladderType));
						}
						/*foreach (var v in u.DownNeighbours) {
							int vertex2=chainStartEnd2Index[new Tuple<char,int,int> (v.SSE.ChainID, v.SSE.Start, v.SSE.End)];
							edges.Add (new Tuple<int,int> (Math.Min(vertex1,vertex2),Math.Max (vertex1,vertex2)));
						}*/
					});
				}
				edges = edges.Distinct().ToList();  // a ladder contained a bulge, there would be listed two edges connecting the strands; this removes the duplicite edges

				SecStrAssignment result = new SecStrAssignment (resultSSEs);
				result.Connectivity = edges;
				result.HBonds = hBonds;
				return result;
			}

			public String GetDescription(){
				return "hydrogen-bond-based method (similar to DSSP, accepted types: E)";
			}

		}

	}
}


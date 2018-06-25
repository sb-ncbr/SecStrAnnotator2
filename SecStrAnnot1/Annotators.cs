using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;

namespace protein
{
	static class Annotators
	{
		public interface IAnnotator{
			List<Tuple<int,int>> GetMatching();
		}

		public class DynProgAnnotator:IAnnotator{
			public AnnotationContext Context{get; private set;}
			public DynProgAnnotator(AnnotationContext context){
				context.ValidateOrdering ();
				this.Context=context;
			}
			public static DynProgAnnotator New(AnnotationContext context){
				return new DynProgAnnotator (context);
			}

			public List<Tuple<int,int>> GetMatching(){
				int m = Context.Templates.Length;
				int n = Context.Candidates.Length;
				double[,] metricMatrix = new double[m, n];
				double[,] dynProgMatrix = new double[m + 1, n + 1];
				int[,] reconstrMatrix = new int[m + 1, n + 1];
				double[] skipTemplatePen = Context.Templates.Select (sse => Context.SkipTemplatePenalty (sse)).ToArray ();
				double[] skipCandidatePen = Context.Candidates.Select (sse => Context.SkipCandidatePenalty (sse)).ToArray ();

				//coding reconstrMatrix: 
				//  -1: optimum was obtained from the left neighbor (skip candidate SSE)
				//   1: optimum was obtained from upper neighbor (skip template SSE)
				//   0: optimum was obtained from diagonal neighbor (pair a template and candidate SSE)

				// Calculation of the matrix of metric.
				for (int i = 0; i < m; i++) {
					for (int j = 0; j < n; j++) {
						metricMatrix [i, j] = Context.MetricToMin (Context.Templates [i], Context.Candidates [j]);
					}
				}

				// Calculation of dynamic programming matrix.
				dynProgMatrix [0, 0] = 0.0;
				for (int i = 1; i < m + 1; i++) {
					reconstrMatrix [i, 0] = 1;
					dynProgMatrix [i, 0] = dynProgMatrix [i - 1, 0] + skipTemplatePen [i - 1];
				}
				for (int j = 1; j < n + 1; j++) {
					reconstrMatrix [0, j] = -1;
					dynProgMatrix [0, j] = dynProgMatrix [0, j - 1] + skipCandidatePen [j - 1];
				}
				for (int i = 1; i < m + 1; i++) {
					for (int j = 1; j < n + 1; j++) {
						SSEInSpace temp = Context.Templates [i - 1];
						SSEInSpace cand = Context.Candidates [j - 1];

						double valueUp = dynProgMatrix [i - 1, j] + skipTemplatePen [i - 1];
						double valueLeft = dynProgMatrix [i, j - 1] + skipCandidatePen [j - 1];

						if (valueUp <= valueLeft) {
							//skip template
							reconstrMatrix [i, j] = 1;
							dynProgMatrix [i, j] = valueUp;
						} else {
							//skip candidate
							reconstrMatrix [i, j] = -1;
							dynProgMatrix [i, j] = valueLeft;
						}

						if (Context.TypeMatching (temp.Type, cand.Type) && (Context.GuideDiscriminator == null || Context.GuideDiscriminator [i - 1, j - 1])) {
							double valueDiag = dynProgMatrix [i - 1, j - 1] + metricMatrix [i - 1, j - 1];
							if (valueDiag <= dynProgMatrix [i, j]) {
								//pair template with candidate
								reconstrMatrix [i, j] = 0;
								dynProgMatrix [i, j] = valueDiag;
							}
						}
					}
				}

				// print dynamic programming matrix and metric matrix
				if (Lib.DoWriteDebug) {
					PrintDynProgAndOrReconstructionMatrix (dynProgMatrix, reconstrMatrix);
					PrintMetricMatrix (Context, metricMatrix, "metric_matrix.tsv");
				}

				// Reconstruction of the best solution.
				List<Tuple<int,int>> matching = new List<Tuple<int, int>> ();
				int row = m;
				int col = n;
				while (row != 0 || col != 0) {
					switch (reconstrMatrix [row, col]) {
					case -1:
						// skipped candidate SSE
						col--;
						break;
					case 1:
						// skipped template SSE
						row--;
						break;
					case 0:
						// paired template and candidate SSE
						matching.Add (new Tuple<int, int> (row - 1, col - 1));
						row--;
						col--;
						break;
					}
				}

				return matching.AsEnumerable ().Reverse ().ToList ();
			}

			private void PrintDynProgAndOrReconstructionMatrix(double[,] dynProgMatrix, int[,] reconstructionMatrix){
				TextWriter w = new StreamWriter (Path.Combine (MainClass.Directory, "dyn_prog_matrix.tsv"));
				w.Write ("\t");
				for (int j = 0; j < Context.Candidates.Length+1; j++) {
					w.Write (j == 0 ? "-\t" : Context.Candidates [j - 1].Label + "\t");
				}
				w.WriteLine ();
				for (int i = 0; i < Context.Templates.Length+1; i++) {
					if (reconstructionMatrix != null) {
						w.Write ("\t");
						for (int j = 0; j < Context.Candidates.Length + 1; j++) {
							w.Write (reconstructionMatrix [i, j] == -1 ? ">" : reconstructionMatrix [i, j] == 1 ? "v" : "\\");
							w.Write ("\t");
						}
						w.WriteLine ();
					}
					if (dynProgMatrix != null) {
						w.Write (i == 0 ? "-\t" : Context.Templates [i - 1].Label + "\t");
						for (int j = 0; j < Context.Candidates.Length + 1; j++) {
							w.Write (dynProgMatrix [i, j].ToString ("0.0") + "\t");
						}
						w.WriteLine ();
					}
				}
				w.Close ();
			}
		}

		public class BranchAndBoundAnnotator:IAnnotator{
			public AnnotationContext Context{get; private set;}
			public BranchAndBoundAnnotator(AnnotationContext context){
				context.ValidateOrdering ();
				context.ValidateBetaGraph ();
				this.Context=context;
			}
			public static BranchAndBoundAnnotator New(AnnotationContext context){
				return new BranchAndBoundAnnotator (context);
			}

			public List<Tuple<int,int>> GetMatching(){
				int m = Context.Templates.Length;
				int n = Context.Candidates.Length;
				double[,] scores = new double[m, n].Fill ((i, j) => 
					Context.TypeMatching (Context.Templates [i].Type, Context.Candidates [j].Type) && (Context.GuideDiscriminator==null || Context.GuideDiscriminator[i,j]) ? 
					Context.SkipTemplatePenalty (Context.Templates [i]) + Context.SkipCandidatePenalty (Context.Candidates [j]) - Context.MetricToMin (Context.Templates [i], Context.Candidates [j]) 
					: 0);
				if (Lib.DoWriteDebug) {
					PrintMetricMatrix (Context, scores, "score_matrix.tsv");
				}
				if (Context.CandidateExclusivity != null)
					PrintCandidateMatrix (Context, Context.CandidateExclusivity, "c_exclusivity.tsv");
				if (Context.CandidateConnectivity != null)
					PrintCandidateMatrix (Context, Context.CandidateConnectivity, "c_connectivity.tsv");

				/*Lib.WriteLineDebug ("Candidate exclusivity:");
				for (int i = 0; i < n; i++) {
					for (int j = i+1; j < n; j++) {
						if (Context.CandidateExclusivity [i, j])
							Lib.WriteLineDebug ("Ex: {0}-{1} / {2}-{3}", Context.Candidates [i].Start, Context.Candidates [i].End, Context.Candidates [j].Start, Context.Candidates [j].End);
					}
				}*/

				DateTime stamp = DateTime.Now;
				Lib.WriteLineDebug ("BranchAndBoundAnnotator.GetAnnotation(): initialized - {0} vs. {1} vertices ({2})", m, n, stamp);

				List<Tuple<int,int>> bestMatching = LibAnnotation.MaxWeightOrderedMatching (m, n, Context.TemplateConnectivity, Context.CandidateConnectivity, scores, Context.TemplateExclusivity,Context.CandidateExclusivity);

				Lib.WriteLineDebug ("BranchAndBoundAnnotator.GetAnnotation(): found matching ({0})", DateTime.Now);
				Lib.WriteLineDebug ("BranchAndBoundAnnotator.GetAnnotation(): time: {0}", DateTime.Now - stamp);


				return bestMatching;
			}
		}

		public class MOMAnnotator:IAnnotator{
			public MOMAnnotationContext MContext{get; private set;}
			public bool SoftOrderConsistency{ get; private set; }

			public MOMAnnotator (AnnotationContext context, bool softOrderConsistency){
				MContext = new MOMAnnotationContext(context);
				SoftOrderConsistency=softOrderConsistency;
			}

			public List<Tuple<int,int>> GetMatching(){
				MContext.Context.ValidateOrdering ();
				int m = MContext.TemplateFST.Length;
				int n = MContext.CandidateFST.Length;

				Func<int,int,double> score = (i, j) =>
					MContext.Context.SkipTemplatePenalty (MContext.Context.Templates [i]) + MContext.Context.SkipCandidatePenalty (MContext.Context.Candidates [j]) - MContext.Context.MetricToMin (MContext.Context.Templates [i], MContext.Context.Candidates [j]);
				Func<int,int,bool> typeMatch = (i, j) =>
					MContext.Context.TypeMatching (MContext.Context.Templates [i].Type, MContext.Context.Candidates [j].Type);
				Func<int,bool> templateIsHelix = (i) =>
					MContext.Context.Templates [MContext.TemplateFST[i].Item1].IsHelix;
				Func<int,int,bool> guideDis = (i, j) =>
					MContext.Context.GuideDiscriminator == null || MContext.Context.GuideDiscriminator [i, j];

				double LADDER_SCORE_SCALE = 1.0;

				double[,] scores = new double[m, n].Fill ((i, j) => 
					MContext.TemplateFST [i].Item3 == MContext.CandidateFST [j].Item3
				                   && typeMatch (MContext.TemplateFST [i].Item1, MContext.CandidateFST [j].Item1)
				                   && typeMatch (MContext.TemplateFST [i].Item2, MContext.CandidateFST [j].Item2)
				                   && guideDis (MContext.TemplateFST [i].Item1, MContext.CandidateFST [j].Item1)
				                   && guideDis (MContext.TemplateFST [i].Item2, MContext.CandidateFST [j].Item2) ?
					(templateIsHelix(i) ? 
						score (MContext.TemplateFST [i].Item1, MContext.CandidateFST [j].Item1) 
						: LADDER_SCORE_SCALE * (score (MContext.TemplateFST [i].Item1, MContext.CandidateFST [j].Item1) + score (MContext.TemplateFST [i].Item2, MContext.CandidateFST [j].Item2))
					)
					: 0);
				if (Lib.DoWriteDebug) {
					PrintMatrix (
						MContext.TemplateFST.Select (t => MContext.Context.Templates [t.Item1].Label + "-" + MContext.Context.Templates [t.Item2].Label).ToArray (),
						MContext.CandidateFST.Select (t => MContext.Context.Candidates [t.Item1].Label + "-" + MContext.Context.Candidates [t.Item2].Label).ToArray (),
						null as object[],
						scores,
						"score_matrix.tsv");
				}
				/*PrintMetricMatrix (Context, scores);
				if (Context.CandidateExclusivity != null)
					PrintCandidateMatrix (Context, Context.CandidateExclusivity, "c_exclusivity.tsv");
				if (Context.CandidateConnectivity != null)
					PrintCandidateMatrix (Context, Context.CandidateConnectivity, "c_connectivity.tsv");*/

				DateTime stamp = DateTime.Now;
				Lib.WriteLineDebug ("{3}.GetAnnotation(): initialized - {0} vs. {1} vertices ({2})",m, n, stamp,this.GetType ().Name);

				List<Tuple<int,int>> bestMOM = LibAnnotation.MaxWeightMixedOrderedMatching (m, n,
					MContext.TemplateFST.Select (t=>t.Item1).ToArray (),
					MContext.TemplateFST.Select (t=>t.Item2).ToArray (),
					MContext.CandidateFST.Select (t=>t.Item1).ToArray (),
					MContext.CandidateFST.Select (t=>t.Item2).ToArray (),
					scores,SoftOrderConsistency).OrderBy (x=>x).ToList ();

				Lib.WriteLineDebug ("MOMAnnotator.GetAnnotation(): found matching ({0})", DateTime.Now);
				Lib.WriteLineDebug ("MOMAnnotator.GetAnnotation(): time: {0}", DateTime.Now - stamp);

				List<Tuple<int,int>> bestMatching = new List<Tuple<int, int>> ();
				foreach (Tuple<int,int> t in bestMOM) {
					Tuple<int,int,int> tLadder = MContext.TemplateFST[t.Item1];
					Tuple<int,int,int> qLadder = MContext.CandidateFST[t.Item2];
					if (tLadder.Item3 == 0) {//helix
						bestMatching.Add (new Tuple<int,int>(tLadder.Item1, qLadder.Item1));
					} else {//ladder
						bestMatching.Add (new Tuple<int,int>(tLadder.Item1, qLadder.Item1));
						bestMatching.Add (new Tuple<int,int>(tLadder.Item2, qLadder.Item2));
					}
				}

				Lib.WriteLineDebug ("Template FST: {0}", MContext.TemplateFST.Select (t=>"("+MContext.Context.Templates[t.Item1].Label+","+MContext.Context.Templates[t.Item2].Label+")").EnumerateWithCommas ());
				Lib.WriteLineDebug ("Candidate FST: {0}", MContext.CandidateFST.Select (t=>"("+MContext.Context.Candidates[t.Item1].Label+","+MContext.Context.Candidates[t.Item2].Label+")").EnumerateWithCommas ());
				Lib.WriteLineDebug ("MOM matching: {0}", bestMOM.EnumerateWithCommas ());
				Lib.WriteLineDebug ("Matching: {0}", bestMatching.Select (t=>"("+MContext.Context.Templates[t.Item1].Label+","+MContext.Context.Candidates[t.Item2].Label+")").EnumerateWithCommas ());

				return bestMatching;
			}

		}

		public class CombinedAnnotator:IAnnotator{
			public AnnotationContext Context{get; private set;}
			public CombinedAnnotator(AnnotationContext context){
				context.ValidateOrdering ();
				context.ValidateBetaGraph ();
				this.Context=context;
			}
			public static CombinedAnnotator New(AnnotationContext context){
				return new CombinedAnnotator (context);
			}

			public List<Tuple<int,int>> GetMatching (){
				DynProgAnnotator dynProgAn = new DynProgAnnotator (Context);
				List<Tuple<int,int>> dynProgMatching = dynProgAn.GetMatching ();
				if (CheckConnectivity (Context, dynProgMatching)) {
					return dynProgMatching;
				} else {
					Lib.WriteLineDebug ("CombinedAnnotator.GetMatching(): Dynamic programming caused connectivity issues. Trying dynamic programming with guide strand matching from branch-and-bound.");
					Lib.Shuffler tShuffler;
					Lib.Shuffler qShuffler;
					IEnumerable<SSEInSpace> filteredTemplates = Context.Templates.WhereAndGetShuffler (sse => sse.IsSheet, out tShuffler);
					IEnumerable<SSEInSpace> filteredCandidates = Context.Candidates.WhereAndGetShuffler (sse => sse.IsSheet, out qShuffler);
					int[,] filteredTemplateConnectivity = tShuffler.ShuffleColumns (tShuffler.ShuffleRows (Context.TemplateConnectivity));
					int[,] filteredCandidateConnectivity = qShuffler.ShuffleColumns (qShuffler.ShuffleRows (Context.CandidateConnectivity));

					AnnotationContext bbContext = new AnnotationContext (Context.MetricToMin, Context.TypeMatching,
						                              Context.SkipTemplatePenalty, Context.SkipCandidatePenalty, 
						                              filteredTemplates, filteredCandidates);
					bbContext.TemplateConnectivity = filteredTemplateConnectivity;
					bbContext.CandidateConnectivity = filteredCandidateConnectivity;

					BranchAndBoundAnnotator bb = new BranchAndBoundAnnotator (bbContext);
					List<Tuple<int,int>> guideMatching = bb.GetMatching ().Select (m => new Tuple<int,int> (tShuffler.OldIndex (m.Item1), qShuffler.OldIndex (m.Item2))).ToList ();
					bool[,] guideDiscriminator = new bool[Context.Templates.Length, Context.Candidates.Length]
						.Fill ((i, j) => !Context.Templates [i].IsSheet && !Context.Candidates [j].IsSheet);
					foreach (var m in guideMatching) {
						guideDiscriminator [m.Item1, m.Item2] = true;
					}
					AnnotationContext dpContext = new AnnotationContext (Context.MetricToMin, Context.TypeMatching,
						                              Context.SkipTemplatePenalty, Context.SkipCandidatePenalty, 
						                              Context.Templates, Context.Candidates);
					dpContext.GuideDiscriminator = guideDiscriminator;
					DynProgAnnotator dp = new DynProgAnnotator (dpContext);
					return dp.GetMatching ();
				}
			}
		}

		public class NiceAnnotatorWrapper{
			const double SUSPICIOUSNESS_THRESHOLD = 1.0;
			const int SUSPICIOUS_JOINING_GAP = 5;

			private IAnnotator inner;
			public AnnotationContext Context{get; private set;}

			public Lib.Shuffler TemplateMapping{ get; private set; }
			public Lib.Shuffler CandidateMapping{ get; private set; }
			private List<Tuple<int,int>> rememberedMatching;
			private SSEInSpace[] rememberedAnnotatedCandidates;
			public bool DoCheckSheetIDConsistency{ get; set; }
			public bool DoRenameSheetIDs{ get; set; }

			public NiceAnnotatorWrapper(AnnotationContext context, Func<AnnotationContext,IAnnotator> innerAnnotatorConstructor){
				Context=context;
				Lib.Shuffler m;
				Context.Templates.OrderAndGetShuffler (out m);
				TemplateMapping=m;
				Context.Candidates.OrderAndGetShuffler(out m);
				CandidateMapping=m;
				DoCheckSheetIDConsistency=true;
				DoRenameSheetIDs=true;

				AnnotationContext innerContext = Context.Copy ();
				innerContext.ApplyShufflers (TemplateMapping,CandidateMapping);
				inner=innerAnnotatorConstructor(innerContext);

				Lib.WriteLineDebug ("Templates: {0}", Context.Templates.Select (s=>s.Start).EnumerateWithCommas ());
				Lib.WriteLineDebug ("Template mapping: {0}", TemplateMapping);
				Lib.WriteLineDebug ("Mapped templates: {0}", innerContext.Templates.Select (s=>s.Start).EnumerateWithCommas ());

			}

			public List<Tuple<int,int>> GetMatching(){
				if (inner == null)
					throw new Exception ("Inner annotator has not been initialized.");
				if (rememberedMatching == null)
					rememberedMatching = inner.GetMatching ().Select (m => new Tuple<int,int> (TemplateMapping.OldIndex (m.Item1), CandidateMapping.OldIndex (m.Item2))).ToList ();
				return rememberedMatching;
			}

			public virtual IEnumerable<SSEInSpace> GetAnnotatedCandidates(){
				if (rememberedAnnotatedCandidates == null) {
					List<Tuple<int,int>> matching = GetMatching ();
					/*Lib.WriteLineDebug ("Matching: {0}",matching.EnumerateWithCommas ());
					Lib.WriteLineDebug ("temps: {0}", matching.Select (t => t.Item1).OrderBy (x => x).EnumerateWithCommas ());
					Lib.WriteLineDebug ("cands: {0}", matching.Select (t => t.Item2).OrderBy (x => x).EnumerateWithCommas ());*/
					if (matching.Select (m => m.Item1).Distinct ().Count () == matching.Count && matching.Select (m => m.Item2).Distinct ().Count () == matching.Count) {
						// each vertex is matched at most once
						Lib.Shuffler annotShuffler = Lib.Shuffler.FromMatching (matching);
						rememberedAnnotatedCandidates = annotShuffler.ShuffleBack (Context.Candidates, () => SSEInSpace.NewNotFound (null))
							.Select ((sse, i) => sse.RelabeledCopy (Context.Templates [i].Label)).ToArray ();
					} else {
						// some vertices are matched more than once
						matching = matching.OrderBy (m => m.Item1).ToList ();
						Dictionary<int,List<int>> multiMatching = new Dictionary<int, List<int>> ();
						foreach (var m in matching)
							multiMatching.MultidictionaryAdd (m.Item1, m.Item2);
						rememberedAnnotatedCandidates = new SSEInSpace[Context.Templates.Length];
						for (int i = 0; i < Context.Templates.Length; i++) {
							if (!multiMatching.ContainsKey (i)) {
								rememberedAnnotatedCandidates [i] = SSEInSpace.NewNotFound (Context.Templates[i].Label);
							} else {
								SSEInSpace[] all = multiMatching [i].Select (j => Context.Candidates [j]).OrderBy (x => x).Distinct ().ToArray ();
								if (all.Length == 1) {
									rememberedAnnotatedCandidates [i] = all [0].RelabeledCopy (Context.Templates [i].Label);
								} else if (all.Length > 1) {
									if (all.Select (sse => sse.ChainID).Distinct ().Count () > 1)
										throw new Exception ("Strands from different chains have been annotated as parts of one strand.");
									for (int j = 0; j < all.Length-1; j++) {
										if(all[j+1].Start-all[j].End-1>=SUSPICIOUS_JOINING_GAP){
											Lib.WriteWarning("Suspicious joining in {0}. Gap between joined SSEs = {1}",Context.Templates [i].Label,all[j+1].Start-all[j].End-1);
										}
									}
									string chainID = all.First ().ChainID;
									int first = all.Select (sse => sse.Start).ArgMin ();
									int last = all.Select (sse => sse.End).ArgMax ();
									rememberedAnnotatedCandidates [i] = new SSEInSpace (new SSE (Context.Templates [i].Label, chainID, all [first].Start, all [last].End, 
										all.Select (sse => sse.Type).Aggregate<char> ((x, y) => MainClass.JoiningTypeCombining (x, y) ?? SSE.NOT_FOUND_TYPE), Context.Templates [i].SheetId), all [first].StartVector, all [last].EndVector);
									rememberedAnnotatedCandidates [i].AddComment ("Created by joining " + all.Count () + " SSEs: " + all.Select (sse => sse.Label).EnumerateWithCommas () + ".");
									foreach (SSEInSpace sse in all) {
										rememberedAnnotatedCandidates [i].AddNestedSSE (sse);
									}
								} else
									throw new Exception ("This should never happen!");
							}
						}
					}
					if (DoCheckSheetIDConsistency || DoRenameSheetIDs) {
						List<Tuple<int,int>> sheetIDMatching;
						CheckSheetIDConsistency (Context, matching, out sheetIDMatching);
						if (DoRenameSheetIDs && sheetIDMatching != null) {
							RenameSheetIDs (rememberedAnnotatedCandidates, sheetIDMatching);
						}
					}
				}
				return rememberedAnnotatedCandidates;
			}

			public virtual List<Tuple<int,int,int>> GetAnnotatedConnectivity(List<Tuple<int,int,int>> connectivity){
				Lib.Shuffler annotShuffler = Lib.Shuffler.FromMatching (GetMatching()).Inverted();
				List<Tuple<int,int,int>> annotConnectivity = annotShuffler.UpdateIndices (connectivity).ToList();
				return annotConnectivity;
			}

			public IEnumerable<T> SelectFromAnnotated<T> (Func<SSEInSpace,SSEInSpace,T> selector){
				return GetAnnotatedCandidates ()
					.Select ((sse, i) => selector (Context.Templates [i], sse));
			}
			public IEnumerable<double> GetMetricList(){
				return SelectFromAnnotated ((t, q) => (q != null && !q.IsNotFound ()) ? Context.MetricToMin (t, q) : 0);
			}

			public void RenameSheetIDs(IEnumerable<SSEInSpace> sses, List<Tuple<int,int>> sheetIDMatching){
				Lib.Shuffler shuffler = Lib.Shuffler.FromMatching (sheetIDMatching);
				foreach (SSEInSpace sse in sses) {
					if (sse.SheetId != null) {
						sse.SheetId = shuffler.OldIndex ((int)sse.SheetId);
					}
				}
			}

			public List<double> GetSuspiciousnessList(){
				List<Tuple<int,int>> matching = GetMatching ();
				double[] result = new double[Context.Templates.Length];
				foreach (var m in matching) {
					SSEInSpace template = Context.Templates [m.Item1];
					SSEInSpace candidate = Context.Candidates [m.Item2];
					double otherMinInRow = Context.Candidates
						.Where ((c, i) => i != m.Item2 && Context.TypeMatching (template.Type, c.Type))
						.Select (c => (double?)Context.MetricToMin (template, c))
						.Min () ?? Double.PositiveInfinity;
					double otherMinInColumn = Context.Templates
						.Where ((t, i) => i != m.Item1 && Context.TypeMatching (t.Type, candidate.Type))
						.Select (t => (double?)Context.MetricToMin (t, candidate))
						.Min () ?? Double.PositiveInfinity;
					double otherMin = Math.Min (otherMinInRow, otherMinInColumn);
					double annotVal = Context.MetricToMin (template, candidate);
					double suspiciousness = annotVal / (annotVal + otherMin);
					result[m.Item1] = suspiciousness;
					if (suspiciousness >= SUSPICIOUSNESS_THRESHOLD) {
						Lib.WriteWarning ("Possibly ambiguous annotation around \"{0}\" (annotated / (alternative+annotated) = {1}).", template.Label, suspiciousness.ToString ("0.00"));
					}
				}
				return result.ToList ();
			}
		}
	
		public class NiceAnnotatorWrapperWithCorrections : NiceAnnotatorWrapper{
			private List<List<String>> corrs; // each list contains [pdbid, label, chainid, startresi, endresi]
			Protein protein;
			public NiceAnnotatorWrapperWithCorrections(
				AnnotationContext context, 
				Func<AnnotationContext,IAnnotator> innerAnnotatorConstructor,
				String correctionsFile,
				String pdbId,
				Protein p)
				: base(context,innerAnnotatorConstructor){
					using (StreamReader r = new StreamReader(correctionsFile)){
						this.corrs = Lib.ReadCSV (r,5,'\t','#').Where (l => l[0]==pdbId).ToList ();
					}
					this.protein = p;
			}
			public override IEnumerable<SSEInSpace> GetAnnotatedCandidates(){
				SSEInSpace[] result = base.GetAnnotatedCandidates ().ToArray ();
				for (int i = 0; i < result.Length; i++) {
					List<String> correction = corrs.Where (l => l [1] == Context.Templates [i].Label).DefaultIfEmpty (null).FirstOrDefault ();
					if (correction != null) {
						SSEInSpace template = Context.Templates [i];
						string chainID = correction [2];
						int start = Int32.Parse (correction [3]);
						int end = Int32.Parse (correction [4]);
						SSE correctedSSE = new SSE (template.Label, chainID, start, end, template.Type, null);
						List<double>[] dump;
						result[i] = (start==0 && end==0) ? 
							SSEInSpace.NewNotFound (template.Label) 
							: LibAnnotation.SSEsAsLineSegments_GeomVersion (protein.GetChain (correctedSSE.ChainID), new SSE[]{ correctedSSE }.ToList (), out dump).First ();
					}
				}
				return result;
			}
		}

		public class AnnotationContext{
			//General settings
			public Func<SSEInSpace,SSEInSpace,double> MetricToMin{ get; private set;}
			public Func<char,char,bool> TypeMatching{ get; private set;}
			public Func<SSEInSpace,double> SkipTemplatePenalty{ get; private set;}
			public Func<SSEInSpace,double> SkipCandidatePenalty{ get; private set;}

			//Templates and candidates
			public SSEInSpace[] Templates{ get; private set;}
			public SSEInSpace[] Candidates{ get; private set;}
			public bool[,] GuideDiscriminator{ get; set; }

			//Beta-graph related context
			public int[,] TemplateConnectivity{ get; set;} //1 = parallel beta-ladder, -1 = antiparallel beta-ladder, 0 = no ladder
			public int[,] CandidateConnectivity{ get; set;} //1 = parallel beta-ladder, -1 = antiparallel beta-ladder, 0 = no ladder
			public bool[,] TemplateExclusivity{ get; set;}
			public bool[,] CandidateExclusivity{ get; set;}

			//Constructors
			public AnnotationContext(
				Func<SSEInSpace,SSEInSpace,double> metricToMin,
				Func<char,char,bool> typeMatching,
				Func<SSEInSpace,double> skipTemplatePenalty,
				Func<SSEInSpace,double> skipCandidatePenalty,
				IEnumerable<SSEInSpace> templates,
				IEnumerable<SSEInSpace> candidates
			){
				this.MetricToMin=metricToMin;
				this.TypeMatching=typeMatching;
				this.SkipTemplatePenalty=skipTemplatePenalty;
				this.SkipCandidatePenalty=skipCandidatePenalty;

				this.Templates=templates.ToArray ();
				this.Candidates=candidates.ToArray ();
			}

			public void InitializeTemplateConnectivity(IEnumerable<Tuple<int,int,int>> templateConnections){
				TemplateConnectivity = new int[Templates.Length, Templates.Length];
				foreach (var connection in templateConnections) {
					TemplateConnectivity [connection.Item1, connection.Item2] = connection.Item3;
					TemplateConnectivity [connection.Item2, connection.Item1] = connection.Item3;
				}
			}

			public void InitializeCandidateConnectivity(IEnumerable<Tuple<int,int,int>> candidateConnections){
				CandidateConnectivity = new int[Candidates.Length, Candidates.Length];
				foreach (var connection in candidateConnections) {
					CandidateConnectivity [connection.Item1, connection.Item2] = connection.Item3;
					CandidateConnectivity [connection.Item2, connection.Item1] = connection.Item3;
				}
			}

			public AnnotationContext Copy(){
				AnnotationContext copy = new AnnotationContext (this.MetricToMin, this.TypeMatching,
					this.SkipTemplatePenalty,this.SkipCandidatePenalty,
					this.Templates,this.Candidates);
				copy.GuideDiscriminator = this.GuideDiscriminator;
				copy.TemplateConnectivity = this.TemplateConnectivity;
				copy.CandidateConnectivity = this.CandidateConnectivity;
				copy.TemplateExclusivity = this.TemplateExclusivity;
				copy.CandidateExclusivity = this.CandidateExclusivity;
				return copy;
			}

			public void ApplyShufflers(Lib.Shuffler templateShuffler, Lib.Shuffler candidateShuffler){
				Templates = templateShuffler.Shuffle (Templates).ToArray ();
				Candidates = candidateShuffler.Shuffle (Candidates).ToArray ();
				if (GuideDiscriminator!=null){
					GuideDiscriminator=templateShuffler.ShuffleRows(candidateShuffler.ShuffleColumns(GuideDiscriminator));
				}

				if (TemplateConnectivity!=null){
					TemplateConnectivity=templateShuffler.ShuffleRowsAndColumns(TemplateConnectivity);
				}

				if (CandidateConnectivity!=null){
					CandidateConnectivity=candidateShuffler.ShuffleRowsAndColumns(CandidateConnectivity);
				}

				if (TemplateExclusivity!=null){
					TemplateExclusivity=templateShuffler.ShuffleRowsAndColumns(TemplateExclusivity);
				}

				if (CandidateExclusivity!=null){
					CandidateExclusivity=candidateShuffler.ShuffleRowsAndColumns(CandidateExclusivity);
				}
			}

			public void ValidateOrdering(){
				if (Templates.Select (x => x.ChainID).Distinct ().Count () > 1)
					Lib.WriteWarning ("{0}: passed template SSEs from more than one chain.", System.Reflection.MethodBase.GetCurrentMethod ().Name);
				if (Candidates.Select (x => x.ChainID).Distinct ().Count () > 1)
					Lib.WriteWarning ("{0}: passed candidate SSEs from more than one chain.", System.Reflection.MethodBase.GetCurrentMethod ().Name);
				
				if (!Templates.IsSorted())
					throw new Exception ("Template SSEs are not ordered.");
				if (!Candidates.IsSorted())
					throw new Exception ("Candidate SSEs are not ordered.");
			}

			public void ValidateBetaGraph(){
				if (TemplateConnectivity==null)
					throw new Exception ("Template connectivity is null.");
				if (CandidateConnectivity==null)
					throw new Exception ("Candidate connectivity is null.");
			}

			public AnnotationContext Ordered(){
				Lib.Shuffler tShuffler;
				Lib.Shuffler cShuffler;
				return Ordered (out tShuffler, out cShuffler);
			}
			public AnnotationContext Ordered(out Lib.Shuffler tShuffler,out Lib.Shuffler cShuffler){
				AnnotationContext result = new AnnotationContext (this.MetricToMin, this.TypeMatching,
					this.SkipTemplatePenalty,this.SkipCandidatePenalty,
					this.Templates.OrderAndGetShuffler (out tShuffler),this.Candidates.OrderAndGetShuffler (out cShuffler));
				result.GuideDiscriminator = tShuffler.ShuffleRows (cShuffler.ShuffleColumns (this.GuideDiscriminator));
				result.TemplateConnectivity = tShuffler.ShuffleRowsAndColumns (this.TemplateConnectivity);
				result.CandidateConnectivity = cShuffler.ShuffleRowsAndColumns (this.CandidateConnectivity);
				result.TemplateExclusivity = tShuffler.ShuffleRowsAndColumns (this.TemplateExclusivity);
				result.CandidateExclusivity = cShuffler.ShuffleRowsAndColumns (this.CandidateExclusivity);
				return result;
			}

			public AnnotationContext Softened(int maxGap){
				//soft matching - joining candidates 
				ValidateOrdering ();
				ValidateBetaGraph ();

				List<SSEInSpace> jointCandidates = new List<SSEInSpace> ();
				List<Tuple<int,int,int>> jointConnections = new List<Tuple<int, int, int>> ();
				List<Tuple<int,int>> conflicts = new List<Tuple<int, int>> ();
				List<Tuple<int,int>> jointGuideAllowedMatches = new List<Tuple<int, int>>();
				for (int i = 0; i < Candidates.Length-1; i++) {
					SSEInSpace s1 = Candidates [i];
					SSEInSpace s2 = Candidates [i + 1];
					if (s1.ChainID==s2.ChainID && s1.IsSheet && s2.IsSheet && s2.Start - s1.End - 1 <= maxGap) {
						int index = Candidates.Length + jointCandidates.Count;
						SSEInSpace newSSE = new SSEInSpace (new SSE (s1.Label + "+" + s2.Label, s1.ChainID, s1.Start, s2.End, 'E', null), s1.StartVector, s2.EndVector);
						newSSE.AddComment ("Created by joining " + s1.Label + " and " + s2.Label + " with gap " + (s2.Start - s1.End - 1) + ".");
						newSSE.AddNestedSSE (s1);
						newSSE.AddNestedSSE (s2);
						jointCandidates.Add (newSSE);
						if (GuideDiscriminator != null) {
							for (int j = 0; j < Templates.Length; j++) {
								if (GuideDiscriminator [j, i] || GuideDiscriminator [j, i + 1])
									jointGuideAllowedMatches.Add (new Tuple<int, int> (j, index));
							}
						}
						for (int j = 0; j < Candidates.Length; j++) {
							if (CandidateConnectivity [i, j] != 0 && CandidateConnectivity [i + 1, j] != 0 && CandidateConnectivity [i, j] != CandidateConnectivity [i, j])
								throw new Exception ("Joint strand would have both parallel and antiparallel ladder to another strand.");
							if (CandidateConnectivity [i, j] == 1 || CandidateConnectivity [i + 1, j] == 1)
								jointConnections.Add (new Tuple<int, int, int> (index, j, 1));
							if (CandidateConnectivity [i, j] == -1 || CandidateConnectivity [i + 1, j] == -1)
								jointConnections.Add (new Tuple<int, int, int> (index, j, -1));
						}
						conflicts.Add (new Tuple<int, int>(i,index));
						conflicts.Add (new Tuple<int, int>(i+1,index));
					}
				}
				SSEInSpace[] newCandidates = Candidates.Concat (jointCandidates).ToArray ();
				bool[,] newGuideDiscriminator;
				if (GuideDiscriminator != null) {
					newGuideDiscriminator = GuideDiscriminator?.Resized2D (Templates.Length, newCandidates.Length);
					foreach (var allowed in jointGuideAllowedMatches) {
						newGuideDiscriminator [allowed.Item1, allowed.Item2] = true;
					}
				} else {
					newGuideDiscriminator = new bool[Templates.Length, newCandidates.Length].Fill (true);
				}
				int[,] newCandidateConnectivity = CandidateConnectivity.Resized2D (newCandidates.Length, newCandidates.Length);
				foreach (var conn in jointConnections) {
					newCandidateConnectivity [conn.Item1, conn.Item2] = conn.Item3;
					newCandidateConnectivity [conn.Item2, conn.Item1] = conn.Item3;
				}
				bool[,] newCandidateExclusivity = CandidateExclusivity?.Resized2D (newCandidates.Length, newCandidates.Length) ?? new bool[newCandidates.Length, newCandidates.Length];
				foreach (var conf in conflicts) {
					newCandidateExclusivity [conf.Item1, conf.Item2] = true;
					newCandidateExclusivity [conf.Item2, conf.Item1] = true;
				}

				//Extra limitations to soft matching to disable useless joining
				int[] tDegreePar = Enumerable.Range (0, Templates.Length).Select (i => Enumerable.Range (0, Templates.Length).Select (j => TemplateConnectivity [i, j]).Count (x=>x==1)).ToArray ();
				int[] tDegreeAnti = Enumerable.Range (0, Templates.Length).Select (i => Enumerable.Range (0, Templates.Length).Select (j => TemplateConnectivity [i, j]).Count (x=>x==-1)).ToArray ();
				int[] jcDegreePar = Enumerable.Range (Candidates.Length,jointCandidates.Count).Select (i => Enumerable.Range (0, Candidates.Length).Select (j => newCandidateConnectivity [i, j]).Count (x=>x==1)).ToArray ();
				int[] jcDegreeAnti = Enumerable.Range (Candidates.Length,jointCandidates.Count).Select (i => Enumerable.Range (0, Candidates.Length).Select (j => newCandidateConnectivity [i, j]).Count (x=>x==-1)).ToArray ();
				for (int i = 0; i < Templates.Length; i++) {
					for (int j = 0; j < jointCandidates.Count; j++) {
						if (jcDegreePar [j] > tDegreePar [i] || jcDegreeAnti [j] > tDegreeAnti [i])
							newGuideDiscriminator [i, Candidates.Length + j] = false;
					}
				}

				Lib.Shuffler shuffler;
				newCandidates = newCandidates.OrderAndGetShuffler (out shuffler).ToArray ();
				AnnotationContext result = new AnnotationContext (MetricToMin, TypeMatching, SkipTemplatePenalty, SkipCandidatePenalty, Templates, newCandidates);
				result.GuideDiscriminator = newGuideDiscriminator != null ? shuffler.ShuffleColumns (newGuideDiscriminator) : null;
				result.TemplateConnectivity = TemplateConnectivity;
				result.CandidateConnectivity = shuffler.ShuffleRowsAndColumns (newCandidateConnectivity);
				result.TemplateExclusivity = TemplateExclusivity;
				result.CandidateExclusivity = shuffler.ShuffleRowsAndColumns (newCandidateExclusivity);

				return result;
			}

			public AnnotationContext Softened_New(int maxGap){
				//soft matching - joining candidates 
				ValidateOrdering ();
				ValidateBetaGraph ();

				List<Tuple<int,int>> joint = new List<Tuple<int, int>> ();
				//Allow joining only 2 strands
				/*for (int i = 0; i < Candidates.Length - 1; i++) {
					int j = i + 1;
					SSEInSpace s1 = Candidates [i];
					SSEInSpace s2 = Candidates [j];
					if (s1.ChainID == s2.ChainID && s1.IsSheet && s2.IsSheet && s2.Start - s1.End - 1 <= maxGap && CandidateConnectivity[i,j]==0) {
						joint.Add (new Tuple<int, int> (i, j));
					}
				}*/
				//Allow joining more than 2 strands
				for (int i = 0; i < Candidates.Length - 1; i++) {
					for (int j = i + 1; j < Candidates.Length; j++) {
						SSEInSpace s1 = Candidates [i];
						SSEInSpace s2 = Candidates [j];
						if (s1.ChainID == s2.ChainID && s1.IsSheet && s2.IsSheet && s2.Start - s1.End - 1 <= maxGap && Enumerable.Range(i,j-i-1).All (x=>CandidateConnectivity[x,j]==0))
							joint.Add (new Tuple<int, int> (i, j));
						else
							break;
					}
				}

				return SoftenedTemplatesOrCandidates (joint,1).Ordered ();
			}

			public AnnotationContext WithAlternativeTemplates(List<Tuple<String,int,int>> alternatives){
				//soft matching - joining templates 
				ValidateBetaGraph ();

				List<String> labels = alternatives.Select (a => a.Item1).ToList ();
				List<Tuple<int,int>> ranges = alternatives.Select (a => new Tuple<int,int> (a.Item2, a.Item3)).ToList ();

				int mOrig = this.Templates.Length;

				AnnotationContext result = this.SoftenedTemplatesOrCandidates (ranges, 0);
				for (int i = 0; i < labels.Count; i++) {
					result.Templates [mOrig + i] = result.Templates [mOrig + i].RelabeledCopy (labels [i]);
				}
				return result.Ordered ();
			}

			/*, whichGraph==0 --> soften templates, whichGraph==1 --> soften candidates; no ordering apllied! */
			private AnnotationContext SoftenedTemplatesOrCandidates(List<Tuple<int,int>> joinedRanges, int whichGraph){
				bool STRICT = false;
				SSEInSpace[] gSSEs; //unchanged
				SSEInSpace[] hSSEs; //to be softened
				bool[,] ghGuideDiscriminator;
				int[,] gConnectivity;
				int[,] hConnectivity;
				bool[,] gExclusivity;
				bool[,] hExclusivity;

				if (whichGraph == 0) {
					gSSEs = Candidates;
					hSSEs = Templates;
					ghGuideDiscriminator = GuideDiscriminator?.Transposed ();
					gConnectivity = CandidateConnectivity;
					hConnectivity = TemplateConnectivity;
					gExclusivity = CandidateExclusivity;
					hExclusivity = TemplateExclusivity;
				} else if (whichGraph == 1) {
					gSSEs = Templates;
					hSSEs = Candidates;
					ghGuideDiscriminator = GuideDiscriminator;
					gConnectivity = TemplateConnectivity;
					hConnectivity = CandidateConnectivity;
					gExclusivity = TemplateExclusivity;
					hExclusivity=CandidateExclusivity;
				} else {
					throw new Exception ("Invalid argument whichGraph.");
				}

				int m = gSSEs.Length;
				int nOrig = hSSEs.Length;
				int nNew = nOrig + joinedRanges.Count;
				hSSEs = hSSEs.Resized (nNew);
				ghGuideDiscriminator = ghGuideDiscriminator?.Resized2D (m, nNew) ?? new bool[m, nNew].Fill (true);
				hConnectivity = hConnectivity.Resized2D (nNew, nNew);
				hExclusivity = hExclusivity?.Resized2D (nNew, nNew) ?? new bool[nNew, nNew];

				List<Tuple<int,int,int>> jointConnections = new List<Tuple<int, int, int>> ();
				List<Tuple<int,int>> conflicts = new List<Tuple<int, int>> ();
				for (int k = 0; k < joinedRanges.Count; k++) {
					int i = joinedRanges [k].Item1;
					int j = joinedRanges [k].Item2;
					int index = nOrig + k;
					int[] partIndices = Enumerable.Range (i, j - i + 1).ToArray ();
					IEnumerable<SSEInSpace> parts = hSSEs.Skip (i).Take (j - i + 1);
					SSEInSpace newSSE = new SSEInSpace (new SSE (parts.Select (p=>p.Label).EnumerateWithSeparators ("_"), parts.First ().ChainID,parts.First ().Start,parts.Last ().End, 'E', null),
						parts.First ().StartVector, parts.Last ().EndVector);
					newSSE.AddComment ("Created by joining " + parts.EnumerateWithSeparators (" and ") + ".");
					foreach (var p in parts) {
						newSSE.AddNestedSSE (p);
					}
					hSSEs [index] = newSSE;
					if (ghGuideDiscriminator != null) {
						for (int t = 0; t < m; t++) {
							ghGuideDiscriminator [t, index] = partIndices.Any (q => ghGuideDiscriminator [t, q]);
						}
					}
					//Connectivity original-joint
					for (int q = 0; q < nOrig; q++) {
						bool parConn = partIndices.Any (p => hConnectivity [q, p] == 1);
						bool antiConn = partIndices.Any (p => hConnectivity [q, p] == -1);
						if (parConn && antiConn)
							throw new Exception ("Joint strand would have both parallel and antiparallel ladder to another strand.");
						if (parConn || antiConn)
							jointConnections.Add (new Tuple<int, int, int> (index, q, parConn ? 1 : -1));
					}
					conflicts.AddRange (partIndices.Select (p=>new Tuple<int,int> (p,index)));
				}

				//Connectivity joint-joint, exclusivity joint-joint
				for (int k = 0; k < joinedRanges.Count; k++) {
					for (int l = k+1; l < joinedRanges.Count; l++) {
						int kIndex = nOrig + k;
						int lIndex = nOrig + l;
						IEnumerable<int> kPartIndices = Enumerable.Range (joinedRanges[k].Item1,joinedRanges[k].Item2-joinedRanges[k].Item1+1);
						IEnumerable<int> lPartIndices = Enumerable.Range (joinedRanges[l].Item1,joinedRanges[l].Item2-joinedRanges[l].Item1+1);
						bool parConn = kPartIndices.Any (kp=>lPartIndices.Any (lp=>hConnectivity[kp,lp]==1));
						bool antiConn = kPartIndices.Any (kp=>lPartIndices.Any (lp=>hConnectivity[kp,lp]==-1));
						if (parConn && antiConn)
							throw new Exception ("Joint strand would have both parallel and antiparallel ladder to another strand.");
						if (parConn || antiConn)
							jointConnections.Add (new Tuple<int, int, int> (kIndex, lIndex, parConn ? 1 : -1));
						if (joinedRanges [k].Item2 >= joinedRanges [l].Item1)
							conflicts.Add (new Tuple<int,int> (kIndex, lIndex));
					}
				}

				foreach (var conn in jointConnections) {
					hConnectivity [conn.Item1, conn.Item2] = conn.Item3;
					hConnectivity [conn.Item2, conn.Item1] = conn.Item3;
				}
				foreach (var conf in conflicts) {
					hExclusivity [conf.Item1, conf.Item2] = true;
					hExclusivity [conf.Item2, conf.Item1] = true;
				}

				//Extra limitations to soft matching to disable useless joining
				if (STRICT) {
					int[] tDegreePar = Enumerable.Range (0, m).Select (i => Enumerable.Range (0, m).Select (j => gConnectivity [i, j]).Count (x => x == 1)).ToArray ();
					int[] tDegreeAnti = Enumerable.Range (0, m).Select (i => Enumerable.Range (0, m).Select (j => gConnectivity [i, j]).Count (x => x == -1)).ToArray ();
					int[] jcDegreePar = Enumerable.Range (nOrig, joinedRanges.Count).Select (i => Enumerable.Range (0, nOrig).Select (j => hConnectivity [i, j]).Count (x => x == 1)).ToArray ();
					int[] jcDegreeAnti = Enumerable.Range (nOrig, joinedRanges.Count).Select (i => Enumerable.Range (0, nOrig).Select (j => hConnectivity [i, j]).Count (x => x == -1)).ToArray ();
					for (int i = 0; i < m; i++) {
						for (int j = 0; j < joinedRanges.Count; j++) {
							if (jcDegreePar [j] > tDegreePar [i] || jcDegreeAnti [j] > tDegreeAnti [i])
								ghGuideDiscriminator [i, nOrig + j] = false;
						}
					}
				}

				/*Lib.Shuffler shuffler;
				hSSEs = hSSEs.OrderAndGetShuffler (out shuffler).ToArray ();
				ghGuideDiscriminator = shuffler.ShuffleColumns (ghGuideDiscriminator);
				hConnectivity = shuffler.ShuffleRowsAndColumns (hConnectivity);
				hExclusivity = shuffler.ShuffleRowsAndColumns (hExclusivity);*/

				if (whichGraph == 0) {
					AnnotationContext result = new AnnotationContext (MetricToMin, TypeMatching, SkipTemplatePenalty, SkipCandidatePenalty, hSSEs, gSSEs);
					result.GuideDiscriminator = ghGuideDiscriminator.Transposed ();
					result.TemplateConnectivity = hConnectivity;
					result.CandidateConnectivity = gConnectivity;
					result.TemplateExclusivity = hExclusivity;
					result.CandidateExclusivity = gExclusivity;
					return result;
				} else if (whichGraph == 1) {
					AnnotationContext result = new AnnotationContext (MetricToMin, TypeMatching, SkipTemplatePenalty, SkipCandidatePenalty, gSSEs, hSSEs);
					result.GuideDiscriminator = ghGuideDiscriminator;
					result.TemplateConnectivity = gConnectivity;
					result.CandidateConnectivity = hConnectivity;
					result.TemplateExclusivity = gExclusivity;
					result.CandidateExclusivity = hExclusivity;
					return result;
				} else {
					throw new Exception ("Invalid argument whichGraph.");
				}
			}
		}

		public class MOMAnnotationContext{
			public AnnotationContext Context{ get; private set; }
			public Tuple<int,int,int>[] TemplateFST{ get; private set; } // tuple<First,Second,Type>
			public Tuple<int,int,int>[] CandidateFST{ get; private set; }
			public MOMAnnotationContext(AnnotationContext inner){
				inner.ValidateOrdering ();
				Context=inner;
				List<Tuple<int,int,int>> templateConnections=new List<Tuple<int, int, int>>();
				for (int i = 0; i < Context.Templates.Length; i++) {
					for (int j = i+1; j < Context.Templates.Length; j++) {
						if (Context.TemplateConnectivity[i,j]!=0){
							templateConnections.Add (new Tuple<int, int, int>(i,j,Context.TemplateConnectivity[i,j]));
						}
					}
				}
				List<Tuple<int,int,int>> candidateConnections=new List<Tuple<int, int, int>>();
				for (int i = 0; i < Context.Candidates.Length; i++) {
					for (int j = i+1; j < Context.Candidates.Length; j++) {
						if (Context.CandidateConnectivity[i,j]!=0){
							candidateConnections.Add (new Tuple<int, int, int>(i,j,Context.CandidateConnectivity[i,j]));
						}
					}
				}

				TemplateFST = Context.Templates.IndicesWhere (sse=>sse.IsHelix).Select (i=>new Tuple<int,int,int> (i,i,0)) //helix nodes
					.Concat (templateConnections.Select (t=>new Tuple<int,int,int>(Math.Min(t.Item1,t.Item2),Math.Max(t.Item1,t.Item2),t.Item3))) //ladder nodes
					.ToArray ();
				CandidateFST = Context.Candidates.IndicesWhere (sse=>sse.IsHelix).Select (i=>new Tuple<int,int,int> (i,i,0)) //helix nodes
					.Concat (candidateConnections.Select (t=>new Tuple<int,int,int>(Math.Min(t.Item1,t.Item2),Math.Max(t.Item1,t.Item2),t.Item3))) //ladder nodes
					.ToArray ();
			}
		}

		public static bool CheckConnectivity (AnnotationContext context, List<Tuple<int,int>> matching){
			return matching.All (m1 => matching.All (m2 => context.TemplateConnectivity [m1.Item1, m2.Item1] == context.CandidateConnectivity [m1.Item2, m2.Item2]));
		}

		/**Return true if sheet IDs are consistent.*/
		public static bool CheckSheetIDConsistency(AnnotationContext context, List<Tuple<int,int>> matching, out List<Tuple<int,int>> sheetIDMatching){
			if (context.Templates.Any (sse => sse.IsSheet && sse.SheetId == null) 
				|| context.Candidates.Any (sse => sse.IsSheet && sse.SheetId == null)) {
				//skip checking (some strand miss sheet ID)
				Lib.WriteWarning ("Cannot check sheet ID consistency because some strand do not have sheet ID.");
				sheetIDMatching = null;
				return true;
			} else {
				Dictionary<int?,List<int?>> dictSheetIdByTemplate = new Dictionary<int?,List<int?>> ();
				Dictionary<int?,List<int?>> dictSheetIdByCandidate = new Dictionary<int?,List<int?>> ();
				foreach (var m in matching) {
					SSE template = context.Templates [m.Item1];
					SSE candidate = context.Candidates [m.Item2];
					if (template.IsSheet && candidate.IsSheet) {
						dictSheetIdByTemplate.MultidictionaryAdd (template.SheetId, candidate.SheetId);
						dictSheetIdByCandidate.MultidictionaryAdd (candidate.SheetId, template.SheetId);
					}
				}
				if (dictSheetIdByTemplate.Values.Any (l => l.Distinct ().Count () > 1) 
					|| dictSheetIdByCandidate.Values.Any (l => l.Distinct ().Count () > 1)) {
					// sheet ID mismatch found
					Lib.WriteWarning ("Some beta-strands caused a sheet ID mismatch!");
					foreach (var kv in dictSheetIdByTemplate.Where (kv => kv.Value.Distinct ().Count () > 1)) {
						Lib.WriteWarning ("Template sheet {0} matched to query sheets {1}. ", kv.Key, kv.Value.Distinct ().EnumerateWithCommas ());
					}
					foreach (var kv in dictSheetIdByCandidate.Where (kv => kv.Value.Distinct ().Count () > 1)) {
						Lib.WriteWarning ("Query sheet {0} matched to template sheets {1}. ", kv.Key, kv.Value.Distinct ().EnumerateWithCommas ());
					}
					sheetIDMatching = null;
					return false;
				} else {
					sheetIDMatching = dictSheetIdByTemplate.Select (kv => new Tuple<int,int> ((int)kv.Key, (int)kv.Value [0])).ToList ();
					return true;
				}
			}
		}

		private static void PrintMetricMatrix<T>(AnnotationContext context, T[,] matrix, String filename){
			PrintMatrix (context.Templates.Select (sse => sse.Label).ToArray (), 
				context.Candidates.Select (sse => sse.Label).ToArray (), 
				context.Templates.Select (sse => (sse.EndVector - sse.StartVector).Size).ToArray (),
				matrix, filename);
		}

		private static void PrintTemplateMatrix<T>(AnnotationContext context, T[,] matrix, String filename){
			String[] names = context.Templates.Select (sse => sse.Label).ToArray ();
			PrintMatrix (names, names, null as object[], matrix, filename);
		}

		private static void PrintCandidateMatrix<T>(AnnotationContext context, T[,] matrix, String filename){
			String[] names = context.Candidates.Select (sse => sse.Label).ToArray ();
			PrintMatrix (names, names, null as object[], matrix, filename);
		}

		public static void PrintCrossMatrix<T>(AnnotationContext context, T[,] matrix, String filename){
			String[] tNames = context.Templates.Select (sse => sse.Label).ToArray ();
			String[] cNames = context.Candidates.Select (sse => sse.Label).ToArray ();
			PrintMatrix (tNames, cNames, null as object[], matrix, filename);
		}

		private static void PrintMatrix<T,T2>(String[] rowNames, String[] colNames, T2[] rowExtras, T[,] matrix, String filename){
			TextWriter w = new StreamWriter (Path.Combine (MainClass.Directory, filename));
			int m = rowNames.Length;
			int n = colNames.Length;
			bool extras = rowExtras != null;
			if (matrix.GetLength (0) != m || matrix.GetLength (1) != n || (extras&& rowExtras.Length!=m))
				throw new Exception ("Wrong dimensions.");
			w.Write ("\t");
			if (extras)
				w.Write ("<extra_info>\t");
			for (int j = 0; j < n; j++) {
				w.Write (colNames [j] + "\t");
			}
			w.WriteLine ();
			for (int i = 0; i < m; i++) {
				w.Write (rowNames [i] + "\t");
				if (extras)
					w.Write (rowExtras[i].ToString () + "\t");
				for (int j = 0; j < n; j++) {
					w.Write (matrix[i,j].ToString (/*"0.0"*/) + "\t");
				}
				w.WriteLine ();
			}
			w.WriteLine ();
			w.WriteLine ();
			w.Close ();
		}
	}
}


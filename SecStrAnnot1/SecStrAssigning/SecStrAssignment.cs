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
	public class SecStrAssignment {
		public List<Sse> SSEs { get; set; }
		public List<(int strand0, int strand1, int type)> Connectivity { get; set; } //each element is a triple <strand1,strand2,ladderType>, where ladderType = 1 for parallel, ladderType = -1 for antiparallel
		public List<(Residue donor, Residue acceptor)> HBonds { get; set; }
		public List<(String mergedLabel, int first, int last)> MergeableSSEs { get; set; } //each element is a triple <label,sse1,sse2> where sse1/sse2 are indices of the first/last SSE to be merged and label is the label of resulting merged SSE

		public SecStrAssignment(IEnumerable<Sse> sses){
			this.SSEs=sses.ToList ();
			this.Connectivity = new List<(int, int, int)>();
			this.HBonds = new List<(Residue, Residue)>();
			this.MergeableSSEs = new List<(String, int, int)>();
		}

		public static SecStrAssignment Combine(SecStrAssignment ass1, SecStrAssignment ass2){
			Lib.Shuffler shuffler;
			SecStrAssignment result = new SecStrAssignment (ass1.SSEs.ConcatAndGetShuffler (ass2.SSEs, out shuffler).ToList ());
			result.Connectivity = ass1.Connectivity.Concat (shuffler.UpdateIndices (ass2.Connectivity)).ToList ();
			result.HBonds = ass1.HBonds.Concat (ass2.HBonds).ToList ();
			return result;
		}
		public static SecStrAssignment Filter(SecStrAssignment ass, Func<Sse,bool> predicate){
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
}


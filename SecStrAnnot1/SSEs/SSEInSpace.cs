using System;
using System.Collections.Generic;
using System.Linq;

using protein.Libraries;
using protein.Geometry;

namespace protein.SSEs
{
	public class SSEInSpace : SSE
	{
		public Point StartVector { get; private set; } // TODO rename to StartPoint
		public Point EndVector { get; private set; }

		public SSEInSpace (SSE sse, Point startVector, Point endVector) 
			: base(sse)
		{
			StartVector = startVector;
			EndVector = endVector;
		}

		public SSEInSpace (SSE sse, (Point, Point) startEndVectors) 
			: this(sse, startEndVectors.Item1, startEndVectors.Item2){ }

		public static new SSEInSpace NewNotFound(String label){
			SSE notFoundSSE= new SSE(label,NOT_FOUND_CHAIN,NOT_FOUND_START, NOT_FOUND_END,NOT_FOUND_TYPE,null);
			notFoundSSE.AddComment ("Not found.");
			return new SSEInSpace (notFoundSSE, new Point(Vector.ZERO), new Point(Vector.ZERO));
		}

		public LineSegment LineSegment => new LineSegment(StartVector, EndVector);

		/** Return an identical SSEInSpace, but with a different label. */
		public new SSEInSpace RelabeledCopy (String newLabel)
		{
			return new SSEInSpace (base.RelabeledCopy(newLabel), StartVector, EndVector);
		}
		public new SSEInSpace RelabeledCopy (String newLabel, String newColor)
		{
			return new SSEInSpace (base.RelabeledCopy(newLabel, newColor), StartVector, EndVector);
		}

		// public static SSEInSpace Join(SSEInSpace first, SSEInSpace second, String comment){
		// 	return new SSEInSpace (SSE.Join(first,second,comment), first.StartVector, second.EndVector);
		// }

		public static SSEInSpace Join(params SSEInSpace[] sses){
			SSEInSpace first = sses[sses.Select(sse => sse.Start).ArgMin()];
			SSEInSpace last = sses[sses.Select(sse => sse.End).ArgMax()];
			SSE newSSE = SSE.Join(sses);
			SSEInSpace newSSEInSpace = new SSEInSpace(newSSE, first.StartVector, last.EndVector);
			return newSSEInSpace;
		}
	}
}


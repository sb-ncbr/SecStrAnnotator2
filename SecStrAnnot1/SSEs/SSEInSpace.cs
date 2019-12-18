using System;
using System.Collections.Generic;
using System.Linq;

using protein.Libraries;
using protein.Geometry;

namespace protein.SSEs
{
	public class SSEInSpace : SSE
	{
		public Point StartPoint { get; private set; }
		public Point EndPoint { get; private set; }

		public SSEInSpace (SSE sse, Point startPoint, Point endPoint) 
			: base(sse)
		{
			StartPoint = startPoint;
			EndPoint = endPoint;
		}

		public SSEInSpace (SSE sse, (Point, Point) startEndPoints) 
			: this(sse, startEndPoints.Item1, startEndPoints.Item2){ }

		public static new SSEInSpace NewNotFound(String label){
			SSE notFoundSSE= new SSE(label,NOT_FOUND_CHAIN,NOT_FOUND_START, NOT_FOUND_END,NOT_FOUND_TYPE,null);
			notFoundSSE.AddComment ("Not found.");
			return new SSEInSpace (notFoundSSE, new Point(Vector.ZERO), new Point(Vector.ZERO));
		}

		public LineSegment LineSegment => new LineSegment(StartPoint, EndPoint);

		/** Return an identical SSEInSpace, but with a different label. */
		public new SSEInSpace RelabeledCopy (String newLabel)
		{
			return new SSEInSpace (base.RelabeledCopy(newLabel), StartPoint, EndPoint);
		}
		public new SSEInSpace RelabeledCopy (String newLabel, String newColor)
		{
			return new SSEInSpace (base.RelabeledCopy(newLabel, newColor), StartPoint, EndPoint);
		}

		// public static SSEInSpace Join(SSEInSpace first, SSEInSpace second, String comment){
		// 	return new SSEInSpace (SSE.Join(first,second,comment), first.StartVector, second.EndVector);
		// }

		public static SSEInSpace Join(params SSEInSpace[] sses){
			SSEInSpace first = sses[sses.Select(sse => sse.Start).ArgMin()];
			SSEInSpace last = sses[sses.Select(sse => sse.End).ArgMax()];
			SSE newSSE = SSE.Join(sses);
			SSEInSpace newSSEInSpace = new SSEInSpace(newSSE, first.StartPoint, last.EndPoint);
			return newSSEInSpace;
		}
	}
}


using System;
using System.Collections.Generic;
using System.Linq;

using protein.Libraries;
using protein.Geometry;

namespace protein.Sses {
    public class SseInSpace : Sse {
        public Point StartPoint { get; private set; }
        public Point EndPoint { get; private set; }
        public Vector? MinorAxis { get; set; }

        public SseInSpace(Sse sse, Point startPoint, Point endPoint)
            : base(sse) {
            StartPoint = startPoint;
            EndPoint = endPoint;
        }

        public SseInSpace(Sse sse, (Point, Point) startEndPoints)
            : this(sse, startEndPoints.Item1, startEndPoints.Item2) { }

        public static new SseInSpace NewNotFound(String label) {
            Sse notFoundSSE = new Sse(label, NOT_FOUND_CHAIN, NOT_FOUND_START, NOT_FOUND_END, SseType.NOT_FOUND_TYPE, null);
            notFoundSSE.AddComment("Not found.");
            return new SseInSpace(notFoundSSE, new Point(Vector.ZERO), new Point(Vector.ZERO));
        }

        public LineSegment LineSegment => new LineSegment(StartPoint, EndPoint);

        /** Return an identical SSEInSpace, but with a different label. */
        public new SseInSpace RelabeledCopy(String newLabel) {
            SseInSpace copy = new SseInSpace(base.RelabeledCopy(newLabel), StartPoint, EndPoint);
            copy.MinorAxis = MinorAxis;
            return copy;
        }
        public new SseInSpace RelabeledCopy(String newLabel, String newColor, String newRainbow) {
            SseInSpace copy = new SseInSpace(base.RelabeledCopy(newLabel, newColor, newRainbow), StartPoint, EndPoint);
            copy.MinorAxis = MinorAxis;
            return copy;
        }

        public static SseInSpace Join(params SseInSpace[] sses) {
            SseInSpace first = sses[sses.Select(sse => sse.Start).ArgMin()];
            SseInSpace last = sses[sses.Select(sse => sse.End).ArgMax()];
            Sse newSSE = Sse.Join(sses);
            SseInSpace newSSEInSpace = new SseInSpace(newSSE, first.StartPoint, last.EndPoint);
            if (sses.All(sse => sse.MinorAxis.HasValue)) {
                newSSEInSpace.MinorAxis = LibAnnotation.CombineMinorAxes(sses.Select(sse => sse.MinorAxis.Value).ToList(), last.EndPoint - first.StartPoint, null);
            }
            return newSSEInSpace;
        }
    }
}


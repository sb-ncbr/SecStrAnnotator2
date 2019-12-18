using System;

namespace protein.Geometry
{
    /** Represents a line segment in 3D given by its start and end point. */
    public class LineSegment
    {
        public Vector Start { get; set; }
        public Vector End { get; set; }
        public LineSegment(Vector start, Vector end)
        {
            Start = start;
            End = end;
        }
        public LineSegment(Point start, Point end) : this(start.Vector, end.Vector) { }
        public Vector Direction { get { return (End - Start).Normalize(); } }
        /** Returns a line that is created by enlongation of this line segment. */
        public Line ToLine()
        {
            return new Line(Start, Direction);
        }
    }
}


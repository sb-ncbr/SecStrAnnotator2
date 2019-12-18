using System;

namespace protein.Geometry
{
    /** Represents a line segment in 3D given by its start and end point. */
    public struct LineSegment
    {
        public readonly Point Start;
        public readonly Point End;

        public LineSegment(Point start, Point end)
        {
            Start = start;
            End = end;
        }

        public Vector Direction => (End - Start).Normalize();
        
        /** Returns a line that is created by enlongation of this line segment. */
        public Line ToLine() => new Line(Start, Direction);

    }
}


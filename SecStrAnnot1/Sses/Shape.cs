using System;
using System.Collections.Generic;
using System.Linq;

namespace protein.Sses
{
	public class Shape{
		public Matrix Points { get; private set; }
		public Matrix Origin { get; private set; }
		public Matrix Axis { get; private set; }
		public Matrix SecondaryAxis { get; private set; }
		public Shape (Matrix points, Matrix origin, Matrix axis, Matrix secondaryAxis)
		{
			this.Points=points;
			this.Origin=origin;
			this.Axis=axis;
			this.SecondaryAxis=secondaryAxis;
		}
	}
}


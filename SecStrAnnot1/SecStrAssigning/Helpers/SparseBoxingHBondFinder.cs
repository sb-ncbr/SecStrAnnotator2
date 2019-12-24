using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;

using Cif.Components;
using protein.Geometry;

namespace protein.SecStrAssigning.Helpers
{
    class SparseBoxingHBondFinder : SimpleHBondFinder
    {
        private const double BOX_SIZE = 10;
        private Dictionary<(int x, int y, int z), List<int>> boxes;
        private Dictionary<(int x, int y, int z), List<int>> neighborBoxes;
        private Point corner;

        public SparseBoxingHBondFinder(IEnumerable<Residue> residues_, double dsspEnergyCutoff)
            : base(residues_, dsspEnergyCutoff)
        {
            boxes = new Dictionary<(int, int, int), List<int>>();
            neighborBoxes = new Dictionary<(int, int, int), List<int>>();

            if (residues.Length > 0)
            {
                corner = new Point(vecCA.Select(v => v.X).Min(), vecCA.Select(v => v.Y).Min(), vecCA.Select(v => v.Z).Min());
                Point otherCorner = new Point(vecCA.Select(v => v.X).Max(), vecCA.Select(v => v.Y).Max(), vecCA.Select(v => v.Z).Max());
                Vector relativeSize = (otherCorner - corner) / BOX_SIZE;
                for (int i = 0; i < vecCA.Length; i++)
                {
                    (int, int, int) indices = GetBoxIndices(i);
                    if (!boxes.ContainsKey(indices))
                        boxes[indices] = new List<int>();
                    boxes[indices].Add(i);
                }
                foreach (var indices in boxes.Keys)
                {
                    neighborBoxes[indices] = MergeNeighborBoxes(indices.x, indices.y, indices.z);
                }
            }
            else
            {
                corner = new Point(Vector.ZERO);
            }
        }

        private (int xIndex, int yIndex, int zIndex) GetBoxIndices(int res)
        {
            Vector relativePosition = (vecCA[res] - corner) / BOX_SIZE;
            return ((int)Math.Floor(relativePosition.X), (int)Math.Floor(relativePosition.Y), (int)Math.Floor(relativePosition.Z));
        }

        private List<int> GetNeighborBox(int res)
        {
            (int xIndex, int yIndex, int zIndex) = GetBoxIndices(res);
            return neighborBoxes[(xIndex, yIndex, zIndex)];
        }

        private List<int> MergeNeighborBoxes(int xIndex, int yIndex, int zIndex)
        {
            List<List<int>> resultBoxes = new List<List<int>>();
            for (int i = xIndex - 1; i <= xIndex + 1; i++)
            {
                for (int j = yIndex - 1; j <= yIndex + 1; j++)
                {
                    for (int k = zIndex - 1; k <= zIndex + 1; k++)
                    {
                        if (boxes.ContainsKey((i, j, k)))
                        {
                            resultBoxes.Add(boxes[(i, j, k)]);
                        }
                    }
                }
            }
            return resultBoxes.SelectMany(b => b).ToList();
        }

        public override List<int> FindHAcceptors(int res)
        {
            return GetNeighborBox(res).Where(s => IsHBond(res, s)).ToList();
        }

        public override List<int> FindHDonors(int res)
        {
            return GetNeighborBox(res).Where(s => IsHBond(s, res)).ToList();
        }

    }
}


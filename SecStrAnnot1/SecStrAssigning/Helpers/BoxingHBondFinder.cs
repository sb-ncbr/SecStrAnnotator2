using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;

using Cif.Components;
using protein.Geometry;

namespace protein.SecStrAssigning.Helpers
{
    class BoxingHBondFinder : SimpleHBondFinder
    {
        private const double BOX_SIZE = 10;
        private List<int>[,,] boxes;
        private List<int>[,,] neighborBoxes;
        private Point corner;
        private readonly (int inX, int inY, int inZ) nBoxes;

        public BoxingHBondFinder(IEnumerable<Residue> residues_, double dsspEnergyCutoff)
            : base(residues_, dsspEnergyCutoff)
        {
            corner = new Point(vecCA.Select(v => v.X).Min(), vecCA.Select(v => v.Y).Min(), vecCA.Select(v => v.Z).Min());
            Point otherCorner = new Point(vecCA.Select(v => v.X).Max(), vecCA.Select(v => v.Y).Max(), vecCA.Select(v => v.Z).Max());
            Vector relativeSize = (otherCorner - corner) / BOX_SIZE;
            nBoxes = ((int)Math.Floor(relativeSize.X) + 1, (int)Math.Floor(relativeSize.Y) + 1, (int)Math.Floor(relativeSize.Z) + 1);
            boxes = new List<int>[nBoxes.inX, nBoxes.inY, nBoxes.inZ];
            for (int i = 0; i < nBoxes.inX; i++)
            {
                for (int j = 0; j < nBoxes.inY; j++)
                {
                    for (int k = 0; k < nBoxes.inZ; k++)
                    {
                        boxes[i, j, k] = new List<int>();
                    }
                }
            }
            for (int i = 0; i < vecCA.Length; i++)
            {
                GetBox(i).Add(i);
            }
            neighborBoxes = new List<int>[nBoxes.inX, nBoxes.inY, nBoxes.inZ];
            for (int i = 0; i < nBoxes.inX; i++)
            {
                for (int j = 0; j < nBoxes.inY; j++)
                {
                    for (int k = 0; k < nBoxes.inZ; k++)
                    {
                        neighborBoxes[i, j, k] = MergeNeighborBoxes(i, j, k);
                    }
                }
            }
            Console.WriteLine("NumBoxes {0}", nBoxes.inX * nBoxes.inY * nBoxes.inZ);
        }
        private (int xIndex, int yIndex, int zIndex) GetBoxIndices(int res)
        {
            Vector relativePosition = (vecCA[res] - corner) / BOX_SIZE;
            return ((int)Math.Floor(relativePosition.X), (int)Math.Floor(relativePosition.Y), (int)Math.Floor(relativePosition.Z));
        }
        private List<int> GetBox(int res)
        {
            (int xIndex, int yIndex, int zIndex) = GetBoxIndices(res);
            return boxes[xIndex, yIndex, zIndex];
        }
        private List<int> GetNeighborBox(int res)
        {
            (int xIndex, int yIndex, int zIndex) = GetBoxIndices(res);
            return neighborBoxes[xIndex, yIndex, zIndex];
        }
        private List<int> MergeNeighborBoxes(int xIndex, int yIndex, int zIndex)
        {
            List<List<int>> result = new List<List<int>>();
            for (int i = Math.Max(0, xIndex - 1); i <= Math.Min(nBoxes.inX - 1, xIndex + 1); i++)
            {
                for (int j = Math.Max(0, yIndex - 1); j <= Math.Min(nBoxes.inY - 1, yIndex + 1); j++)
                {
                    for (int k = Math.Max(0, zIndex - 1); k <= Math.Min(nBoxes.inZ - 1, zIndex + 1); k++)
                    {
                        result.Add(boxes[i, j, k]);
                    }
                }
            }
            return result.SelectMany(b => b).ToList();
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


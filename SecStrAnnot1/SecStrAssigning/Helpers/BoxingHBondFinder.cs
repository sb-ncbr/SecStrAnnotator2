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
        private int[] nBoxes; // TODO change to ValueTuple!
        public BoxingHBondFinder(IEnumerable<Residue> residues_, double dsspEnergyCutoff)
            : base(residues_, dsspEnergyCutoff)
        {
            corner = new Point(vecCA.Select(v => v.X).Min(), vecCA.Select(v => v.Y).Min(), vecCA.Select(v => v.Z).Min());
            Point otherCorner = new Point(vecCA.Select(v => v.X).Max(), vecCA.Select(v => v.Y).Max(), vecCA.Select(v => v.Z).Max());
            Vector relativeSize = (otherCorner - corner) / BOX_SIZE;
            nBoxes = new int[] { (int)Math.Floor(relativeSize.X) + 1, (int)Math.Floor(relativeSize.Y) + 1, (int)Math.Floor(relativeSize.Z) + 1 };
            // nBoxes = new int[3];
            // for (int i = 0; i < 3; i++)
            // {
            //     nBoxes[i] = (int)Math.Floor((otherCorner[i] - corner[i]) / BOX_SIZE) + 1;
            // }
            boxes = new List<int>[nBoxes[0], nBoxes[1], nBoxes[2]];
            for (int i = 0; i < nBoxes[0]; i++)
            {
                for (int j = 0; j < nBoxes[1]; j++)
                {
                    for (int k = 0; k < nBoxes[2]; k++)
                    {
                        boxes[i, j, k] = new List<int>();
                    }
                }
            }
            for (int i = 0; i < vecCA.Length; i++)
            {
                GetBox(i).Add(i);
            }
            neighborBoxes = new List<int>[nBoxes[0], nBoxes[1], nBoxes[2]];
            for (int i = 0; i < nBoxes[0]; i++)
            {
                for (int j = 0; j < nBoxes[1]; j++)
                {
                    for (int k = 0; k < nBoxes[2]; k++)
                    {
                        neighborBoxes[i, j, k] = MergeNeighborBoxes(new int[] { i, j, k });
                    }
                }
            }
        }
        private int[] GetBoxIndices(int res)
        {
            Vector relativePosition = (vecCA[res] - corner) / BOX_SIZE;
            return new int[] { (int)Math.Floor(relativePosition.X), (int)Math.Floor(relativePosition.Y), (int)Math.Floor(relativePosition.Z) };

            // int[] indices = new int[3];
            // double[] pos = vecCA[res].AsArray();
            // for (int i = 0; i < 3; i++)
            //     indices[i] = (int)Math.Floor((pos[i] - corner[i]) / BOX_SIZE);
            // return indices;
        }
        private List<int> GetBox(int res)
        {
            int[] indices = GetBoxIndices(res);
            return boxes[indices[0], indices[1], indices[2]];
        }
        private List<int> GetNeighborBox(int res)
        {
            int[] indices = GetBoxIndices(res);
            return neighborBoxes[indices[0], indices[1], indices[2]];
        }
        private List<int> MergeNeighborBoxes(int[] indices)
        {
            List<List<int>> result = new List<List<int>>();
            for (int i = Math.Max(0, indices[0] - 1); i <= Math.Min(nBoxes[0] - 1, indices[0] + 1); i++)
            {
                for (int j = Math.Max(0, indices[1] - 1); j <= Math.Min(nBoxes[1] - 1, indices[1] + 1); j++)
                {
                    for (int k = Math.Max(0, indices[2] - 1); k <= Math.Min(nBoxes[2] - 1, indices[2] + 1); k++)
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


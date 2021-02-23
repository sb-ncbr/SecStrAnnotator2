
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Cif.Components;
using protein.Json;
using protein.Libraries;
using SecStrAnnotator2.Utils;

namespace protein
{
    public static class DumbAligner
    {
        public static Protein DumbAlign(Protein target, Protein mobile, string alignedMobileOutputFile, string alignmentJsonFile = null)
        {
            Atom[] targetAlphas = target.GetResidues().Select(r => r.GetCAlpha()).WhereHasValue().ToArray();
            Atom[] mobileAlphas = mobile.GetResidues().Select(r => r.GetCAlpha()).WhereHasValue().ToArray();
            targetAlphas[0].Position();
            Geometry.Point[] tPositions = targetAlphas.Select(a => a.Position()).ToArray();
            Geometry.Point[] mPositions = mobileAlphas.Select(a => a.Position()).ToArray();
            int targetLength = targetAlphas.Length;
            int mobileLength = mobileAlphas.Length;
            int alignLength = Math.Min(targetLength, mobileLength);
            Matrix trans;
            Matrix rot;
            double rmsd;
            int minShift = Math.Min(0, targetLength - mobileLength);
            int maxShift = Math.Max(0, targetLength - mobileLength);
            List<(double rmsd, int shift, Matrix rot, Matrix trans)> tranformations = new List<(double, int, Matrix, Matrix)>();
            int shift;
            for (shift = minShift; shift <= maxShift; shift++)
            {
                Matrix t;
                Matrix m;
                if (shift >= 0)
                {
                    t = Matrix.FromRows(tPositions[shift..(alignLength + shift)]);
                    m = Matrix.FromRows(mPositions[0..alignLength]);
                }
                else
                {
                    t = Matrix.FromRows(tPositions[0..alignLength]);
                    m = Matrix.FromRows(mPositions[(-shift)..(alignLength - shift)]);
                }
                LibAlgebra.Align(m, t, out trans, out rot, out rmsd, centeredTarget: false);
                tranformations.Add((rmsd, shift, rot, trans));
            }
            (rmsd, shift, rot, trans) = tranformations.Min();
            Protein transformedMobile = mobile.Transform(rot, trans);
            if (alignedMobileOutputFile != null)
            {
                transformedMobile.SaveCif(alignedMobileOutputFile);
            }
            if (alignmentJsonFile != null)
            {
                Residue[] tResidues = target.GetResidues().ToArray();
                Residue[] mResidues = mobile.GetResidues().ToArray();
                var alignedResidues = new List<object>();
                for (int i = 0; i < alignLength; i++)
                {
                    Residue t;
                    Residue m;
                    if (shift >= 0)
                    {
                        t = tResidues[i + shift];
                        m = mResidues[i];
                    }
                    else
                    {
                        t = tResidues[i];
                        m = mResidues[i - shift];
                    }
                    alignedResidues.Add(new List<object> { t.ChainId, t.SeqNumber, m.ChainId, m.SeqNumber });
                }
                JsonValue j = JsonValue.MakeObject();
                j["alignment_length"] = new JsonValue(alignLength);
                j["RMSD"] = new JsonValue(rmsd);
                j["rotation_matrix"] = new JsonValue(PymolTranformationMatrix(rot, trans).Select(x => x as object).ToList());
                j["aligned_residues"] = new JsonValue(alignedResidues);
                using (StreamWriter w = new StreamWriter(alignmentJsonFile))
                {
                    w.Write(j.ToString(Setting.JSON_OUTPUT_MAX_INDENT_LEVEL));
                }
            }
            return transformedMobile;
        }
        private static double[] PymolTranformationMatrix(Matrix R, Matrix t)
        {
            double[] result = new double[]{
                R[0,0], R[1,0], R[2,0], t[0,0],
                R[0,1], R[1,1], R[2,1], t[0,1],
                R[0,2], R[1,2], R[2,2], t[0,2],
                0.0,    0.0,    0.0,    1.0,
            };
            return result;
        }
        
    }
}
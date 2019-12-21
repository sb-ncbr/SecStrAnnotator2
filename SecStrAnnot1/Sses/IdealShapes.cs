using System;
using System.Collections.Generic;
using System.Linq;

namespace protein.Sses
{
    /**This static class provides "ideal" geometries of alpha-helix and beta-sheet obtained from 1OG2 and 5ACM respectively.
	 * These are defined in IdealShapes.InitializeDictionary().*/
    public static class IdealShapes
    {
        private static Dictionary<SseType, Shape> idealShapeDictionary;

        public static Shape GetShape(SseType type)
        {
            if (idealShapeDictionary == null)
            {
                InitializeDictionary();
            }
            return idealShapeDictionary[type];
        }

        private static void InitializeDictionary()
        {
            Shape idealHelixH = new Shape(
                Matrix.CreateByRows(
                    4, 3, new double[] {
                        -1.71821, -1.24365, -2.27248,
                        1.73155, -1.74169, -0.75852,
                        1.71988, 1.74169, 0.76133,
                        -1.73322, 1.24365, 2.26967
                    }),
                Matrix.CreateByRows(
                    1, 3, new double[] { 0.1876871, -0.0014040, 0.0057961 }),
                Matrix.CreateByRows(
                    1, 3, new double[] { 0, 0, 1 }),
                Matrix.CreateByRows(
                    1, 3, new double[] { 1, 0, 0 })
            );

            Shape idealHelixH5 = new Shape(
                Matrix.CreateByRows(
                    5, 3, new double[] {
                        -1.4821238,  -0.8513020,  -3.0297950,
                         1.9676362,  -1.3493420,  -1.5158350,
                         1.9559662,   2.1340380,   0.0040150,
                        -1.4971338,   1.6359980,   1.5123550,
                        -0.9443448,  -1.5693920,   3.0292600 //5th point approximated
					}),
                Matrix.CreateByRows(
                    1, 3, new double[] { 0, 0, 0 }),
                Matrix.CreateByRows(
                    1, 3, new double[] { 0, 0, 1 }),
                Matrix.CreateByRows(
                    1, 3, new double[] { 1, 0, 0 })
            );

            Shape idealHelixH6 = new Shape(
                Matrix.CreateByRows(
                    6, 3, new double[] {
                        -1.18043, 1.96174, -3.78657, //1st point approximated
						-1.71821, -1.24365, -2.27248,
                        1.73155, -1.74169, -0.75852,
                        1.71988, 1.74169, 0.76133,
                        -1.73322, 1.24365, 2.26967,
                        -1.18043, -1.96174, 3.78657 //6th point approximated
					}),
                Matrix.CreateByRows(
                    1, 3, new double[] { 0.1876871, -0.0014040, 0.0057961 }),
                Matrix.CreateByRows(
                    1, 3, new double[] { 0, 0, 1 }),
                Matrix.CreateByRows(
                    1, 3, new double[] { 1, 0, 0 })
            );

            Shape idealSheetE = new Shape(
                Matrix.CreateByRows(
                4, 3, new double[] {
                -0.15076, -0.81533, -5.03094,
                0.22058, 0.88737, -1.68721,
                0.20383, -0.83079, 1.67805,
                -0.27365, 0.75875, 5.04010
            }),
                Matrix.CreateByRows(
                1, 3, new double[] { 4.2728e-02, -7.4701e-03, 8.2424e-04 }),
                Matrix.CreateByRows(
                1, 3, new double[] { 0, 0, 1 }),
                Matrix.CreateByRows(
                1, 3, new double[] { 1, 0, 0 })
            );

            idealShapeDictionary = new Dictionary<SseType, Shape>() { };
            foreach (SseType type in LibSseTypes.ALL_HELIX_TYPES)
            {
                idealShapeDictionary[type] = idealHelixH;
            }
            foreach (SseType type in LibSseTypes.ALL_SHEET_TYPES)
            {
                idealShapeDictionary[type] = idealSheetE;
            }
            // idealShapeDictionary.Add('5', idealHelixH5);
            // idealShapeDictionary.Add('6', idealHelixH6);

            /*idealShapeDictionary = new Dictionary<char, Shape> () {
				{ 'H', idealHelixH },
				{ 'G', idealHelixH },
				{ 'I', idealHelixH },
				{ SSE.MIXED_HELIX_TYPE, idealHelixH },
				{ 'E', idealSheetE }
			};*/
        }
    }
}


using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using Cif;
using Cif.Tables;
using Cif.Filtering;
using Cif.Components;

namespace SecStrAnnotator2
{
    public static class CifWrapperForSecStrAnnot1_New {
        private const string EXCEPTION_MESSAGE = "Cannot fit CIF data to PDB data model: ";
        private const char PDB_DEFAULT_ALT_LOC = ' ';
        private const char PDB_DEFAULT_INS_CODE = ' ';
        private const double DEFAULT_OCCUPANCY = 1.0;
        private const double DEFAULT_TEMP_FACTOR = 0.0;
        private const string DEFAULT_CHARGE = "  ";

        public static Cif.Components.Protein ProteinFromCifFile(string filename) {
            CifPackage package = CifPackage.FromFile(filename);
            if (package.BlockNames.Length < 1) {
                throw new FormatException(EXCEPTION_MESSAGE + "CIF file must contain at least one block");
            }
            ModelCollection models = ModelCollection.FromCifBlock(package.Blocks[0]);
            if (models.Count < 1) {
                throw new FormatException(EXCEPTION_MESSAGE + "CIF category must contain at least one model");
            }
            Model model = models.GetModel(0);
            return ProteinFromCifModel(model);
        }

        public static Cif.Components.Protein ProteinFromCifFile(string filename, string chainId, (int,int)[] resSeqRanges) {
            CifPackage package = CifPackage.FromFile(filename);
            if (package.BlockNames.Length < 1) {
                throw new FormatException(EXCEPTION_MESSAGE + "CIF file must contain at least one block");
            }
            CifCategory category = package.Blocks[0].GetCategory(ModelCollection.CATEGORY_NAME);
            // int[] rows = category.GetItem(ChainTable.ID_COLUMN).GetRowsWith(chainId.ToString()).ToArray(); // filter rows by chain ID
            // int[] seqNumbers = category.GetItem(ResidueTable.SEQ_NUMBER_COLUMN).GetIntegers(rows, 0);
            // rows = rows.Where( (r, i) => resSeqRanges.Any(range => range.Item1 <= seqNumbers[i] && seqNumbers[i] <= range.Item2) ).ToArray(); // filter rows by resi
            Filter filter = Filter.StringEquals(ChainTable.ID_COLUMN, chainId.ToString()) 
                            & Filter.IntegerInRange(ResidueTable.SEQ_NUMBER_COLUMN, resSeqRanges);
            int[] rows = filter.GetFilteredRows(category).ToArray();

            ModelCollection models = ModelCollection.FromCifCategory(category, rows);
            if (models.Count < 1) {
                string rangeString = String.Join(",", resSeqRanges.Select(range => range.Item1 + ":" + range.Item2));
                throw new FormatException(EXCEPTION_MESSAGE + $"Atom selection given by chain ID '{chainId}' and residue ranges '{rangeString}' is empty");
            }
            Model model = models.GetModel(0);
            Protein protein = ProteinFromCifModel(model);
            return protein;
        }

        public static Cif.Components.Protein ProteinFromCifModel(Model model) {
            return new Cif.Components.Protein(model);
        }        
    }
}
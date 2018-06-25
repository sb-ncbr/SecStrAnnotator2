using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using SecStrAnnot2.Cif;
using SecStrAnnot2.Cif.Tables;

namespace SecStrAnnot2
{
    public static class CifWrapperForSecStrAnnot1 {
        private const string EXCEPTION_MESSAGE = "Cannot fit CIF data to PDB data model: ";
        private const char PDB_DEFAULT_ALT_LOC = ' ';
        private const char PDB_DEFAULT_INS_CODE = ' ';
        private const double DEFAULT_OCCUPANCY = 1.0;
        private const double DEFAULT_TEMP_FACTOR = 0.0;
        private const string DEFAULT_CHARGE = "  ";

        public static protein.Protein ProteinFromCifFile(string filename) {
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

        public static protein.Protein ProteinFromCifFile(string filename, string chainId, IEnumerable<Tuple<int,int>> resSeqRanges) {
            CifPackage package = CifPackage.FromFile(filename);
            if (package.BlockNames.Length < 1) {
                throw new FormatException(EXCEPTION_MESSAGE + "CIF file must contain at least one block");
            }
            CifCategory category = package.Blocks[0].GetCategory(ModelCollection.CATEGORY_NAME);
            int[] rows = category.GetItem(ChainTable.ID_COLUMN).GetRowsWith(chainId.ToString()); // filter rows by chain ID
            int[] seqNumbers = category.GetItem(ResidueTable.SEQ_NUMBER_COLUMN).GetIntegers(rows, 0);
            rows = rows.Where( (r, i) => resSeqRanges.Any(range => range.Item1 <= seqNumbers[i] && seqNumbers[i] <= range.Item2) ).ToArray(); // filter rows by resi

            ModelCollection models = ModelCollection.FromCifCategory(category, rows);
            if (models.Count < 1) {
                throw new FormatException(EXCEPTION_MESSAGE + "Atom selection given by chain ID and residue ranges is empty");
            }
            Model model = models.GetModel(0);
            return ProteinFromCifModel(model);
        }

        public static protein.Protein ProteinFromCifModel(Model model) {
            return new protein.Protein(AtomsFromModel(model));
        }

        private static protein.Atom[] AtomsFromModel(Model model) {
            AtomTable atoms = model.Atoms;
            ResidueTable residues = model.Residues;
            ChainTable chains = model.Chains;

            protein.Atom[] atomArray = new protein.Atom[atoms.Count];

            for (int iAtom = 0; iAtom < atoms.Count; iAtom++) {
                int serial;
                try {
                    serial = int.Parse(atoms.Id[iAtom]);
                } catch {
                    throw new FormatException(EXCEPTION_MESSAGE + $"Atom ID must be an integer: {atoms.Id[iAtom]}");
                }
                string elementString = atoms.Element[iAtom];
                string name = atoms.Name[iAtom];
                // if (name.Length < 4) {
                //     if (elementString.Length == 1) {
                //         name = " " + name;
                //     }
                //     name = name.PadRight(4);
                // } else if (name.Length == 4) {
                //     // this is OK, do nothing
                // } else {
                //     throw new FormatException(EXCEPTION_MESSAGE + $"Atom name must have at most 4 characters: {name}");
                // }
                
                string altLocString = atoms.AltLoc[iAtom];
                char altLoc;
                if (altLocString == AtomTable.DEFAULT_ALT_LOC) {
                    altLoc = PDB_DEFAULT_ALT_LOC;
                } else if (altLocString.Length == 1) {
                    altLoc = altLocString[0];
                } else {
                    throw new FormatException(EXCEPTION_MESSAGE + $"Alternative location must be a single character: {altLocString}");
                }
                string resName = residues.Compound[atoms.ResidueIndex[iAtom]];
                string chainIdString = chains.Id[atoms.ChainIndex[iAtom]];
                // char chainId;
                // if (chainIdString.Length == 1) {
                //     chainId = chainIdString[0];
                // } else {
                //     throw new FormatException(EXCEPTION_MESSAGE + $"Chain ID must be a single character: {chainIdString}");
                // }
                int resSeq = residues.SeqNumber[atoms.ResidueIndex[iAtom]];
                char insCode = PDB_DEFAULT_INS_CODE;
                double x = atoms.X[iAtom];
                double y = atoms.Y[iAtom];
                double z = atoms.Z[iAtom];
                double occupancy = DEFAULT_OCCUPANCY;
                double tempFactor = DEFAULT_TEMP_FACTOR;
                // string element;
                // if (elementString.Length <= 2) {
                //     element = elementString.PadLeft(2);
                // } else {
                //     throw new FormatException(EXCEPTION_MESSAGE + $"Element symbol must be have at most 2 characters: {elementString}");
                // }
                string charge = DEFAULT_CHARGE;
                bool isHetatm = atoms.IsHetatm[iAtom];

                atomArray[iAtom] = new protein.Atom(serial, name, altLoc, resName, chainIdString, resSeq, insCode, x, y, z, occupancy, tempFactor, elementString, charge, isHetatm);
            }

            return atomArray;
        }
        
    }
}
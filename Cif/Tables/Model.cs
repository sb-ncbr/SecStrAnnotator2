using System;
using System.Collections.Generic;
using System.Linq;

namespace SecStrAnnot2.Cif.Tables
{
    public class Model
    {
        public int ModelNumber { get; private set; }
        public AtomTable Atoms { get; private set; }
        public ResidueTable Residues { get; private set; }
        public FragmentTable Fragments { get; private set; }
        public ChainTable Chains { get; private set; }
        public EntityTable Entities { get; private set; }

        ///<summary> Not to be called directly! Use ModelCollection.GetModel() or similar.</summary>
        internal Model(CifCategory category, int[] rows, int modelNumber){
            this.ModelNumber = modelNumber;

            // entities
            CifItem entityId = category[EntityTable.KEY_COLUMN];
            int[] atomStartsOfEntities;
            rows = entityId.GetRowsGroupedByValue(rows, out atomStartsOfEntities);

            // chains
            CifItem asymId = category[ChainTable.KEY_COLUMN];
            int[] atomStartsOfChains;
            int[] chainStartsOfEntities;
            rows = asymId.GetRowsGroupedByValueInEachRegion(rows, atomStartsOfEntities, out atomStartsOfChains, out chainStartsOfEntities); //TODO assume grouped

            // fragments + residues
            CifItem residueId = category[ResidueTable.KEY_COLUMN];
            int[] fragmentStartsOfChains;
            int[] residueStartsOfFragments;
            int[] atomStartsOfResidues;
            int[] residueStartsOfChains;
            int[] residueNumbers;
            GetFragmentsAndResidues(residueId, ref rows, atomStartsOfChains, out fragmentStartsOfChains, out residueStartsOfFragments, out atomStartsOfResidues, out residueNumbers);
            rows = residueId.GetRowsGroupedByValueInEachRegion(rows, atomStartsOfChains, out atomStartsOfResidues, out residueStartsOfChains); //TODO replace by sorting by seq_id! + try to assume sorted
            
            // remaining combinations
            int[] residueStartsOfEntities = GetSelectedElements(residueStartsOfChains, chainStartsOfEntities, false);
            int[] fragmentStartsOfEntities = GetSelectedElements(fragmentStartsOfChains, chainStartsOfEntities, false);

            //TODO create Chains, Fragments, Residues, Atoms in the same way + make those tables nested classes of Model and their constructors only available to Model?
            this.Entities = new EntityTable(category, rows, atomStartsOfEntities, residueStartsOfEntities, fragmentStartsOfEntities, chainStartsOfEntities);

            Lib.LogList("Entities.chainStartIndex", Entities.chainStartIndex);
            Lib.LogList("Entities.fragmentStartIndex", Entities.fragmentStartIndex);
            Lib.LogList("Entities.residueStartIndex", Entities.residueStartIndex);
            Lib.LogList("Entities.atomStartIndex", Entities.atomStartIndex);
            // string[] atomIds = category["id"].GetStrings();
            // string[] atomNames = category["label_atom_id"].GetStrings();
            // string[] seqIds = category["label_seq_id"].GetStrings();
            // string[] compIds = category["label_comp_id"].GetStrings();
            // string[] asymIds = category["label_asym_id"].GetStrings();
            // string[] entityIds = category["label_entity_id"].GetStrings();           
            // Console.WriteLine("row \tatom \ta.name\t  comp seq\tasym \tentity");
            // foreach (int iRow in rows){
            //     Console.WriteLine($"{iRow}:\t {atomIds[iRow]}\t {atomNames[iRow]}\t   {compIds[iRow]}  {seqIds[iRow]} \t {asymIds[iRow]}\t {entityIds[iRow]}");
            // }
        }

        private static void GetFragmentsAndResidues(
            CifItem residueNumberItem, 
            ref int[] rows, 
            int[] atomStartsOfChains, 
            out int[] fragmentStartsOfChains,
            out int[] residueStartsOfFragments,
            out int[] atomStartsOfResidues,
            out int[] residueNumbersOfResidues)
        {
            if (rows.Length == 0){
                Lib.WriteWarning("Creating model from 0 atoms.");
            }
            int[] residueNumbers = residueNumberItem.GetIntegers(rows, ResidueTable.DEFAULT_RESIDUE_NUMBER);
            int nChains = atomStartsOfChains.Length - 1; // the last is sentinel
            List<int> fragOfChainList = new List<int>();
            List<int> resOfFragList = new List<int>();
            List<int> atomOfResList = new List<int>();
            int iResidue = 0;
            int iFragment = 0;
            for (int iChain = 0; iChain < nChains; iChain++)
            {
                int startAtom = atomStartsOfChains[iChain];
                int endAtom = atomStartsOfChains[iChain+1];
                // start new residue, new fragment, new chain
                int resNum = residueNumbers[startAtom];
                atomOfResList.Add(startAtom);      
                resOfFragList.Add(iResidue);   
                fragOfChainList.Add(iFragment);
                for (int iAtom = startAtom+1; iAtom < endAtom; iAtom++)
                {
                    int newResNum = residueNumbers[iAtom];
                    if (newResNum == resNum) {
                        // continue residue (do nothing)
                    } else if (newResNum == resNum + 1) {
                        // new residue, continue fragment
                        //TODO start new residue if this is water
                        resNum = newResNum;
                        atomOfResList.Add(iAtom);
                        iResidue++;
                    } else if (newResNum > resNum) {
                        // new residue, new fragment
                        resNum = newResNum;
                        atomOfResList.Add(iAtom);
                        iResidue++;
                        resOfFragList.Add(iResidue);
                        iFragment++;
                    } else {
                        // decreasing residue number => throw exception
                        throw new NotImplementedException("Residue numbers (" + residueNumberItem.FullName + ") are not in increasing order (not supported by current version)");
                    }
                }
                // finish residue, fragment, chain
                iResidue++;
                iFragment++;
            }
            fragmentStartsOfChains = fragOfChainList.Append(iFragment).ToArray();
            residueStartsOfFragments = resOfFragList.Append(iResidue).ToArray();
            atomStartsOfResidues = atomOfResList.Append(atomStartsOfChains[nChains]).ToArray();
            residueNumbersOfResidues = atomOfResList.Select(a => residueNumbers[a]).ToArray();
        }

        internal static int[] GetSelectedElements(int[] elements, int[] indices, bool ignoreLastIndex) {
            int nSelected = ignoreLastIndex ? indices.Length - 1 : indices.Length;
            int[] selected = new int[nSelected];
            for (int i = 0; i < nSelected; i++)
            {
                selected[i] = elements[indices[i]];
            }
            return selected;
        }
    }
}
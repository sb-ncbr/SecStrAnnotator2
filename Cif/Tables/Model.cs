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

            CifItem entityId = category[EntityTable.KEY_COLUMN];
            int[] atomStartsOfEntities;
            rows = entityId.GetRowsGroupedByValue(rows, out atomStartsOfEntities);

            CifItem asymId = category[ChainTable.KEY_COLUMN];
            int[] atomStartsOfChains;
            int[] chainStartsOfEntities;
            rows = asymId.GetRowsGroupedByValueInEachRegion(rows, atomStartsOfEntities, out atomStartsOfChains, out chainStartsOfEntities); //TODO assume grouped

            CifItem residueId = category[ResidueTable.KEY_COLUMN];
            int[] atomStartsOfResidues;
            int[] residueStartsOfChains;
            int[] residueStartsOfEntities;
            rows = residueId.GetRowsGroupedByValueInEachRegion(rows, atomStartsOfChains, out atomStartsOfResidues, out residueStartsOfChains); //TODO replace by sorting by seq_id! + try to assume sorted

            //TODO fragments

            string[] atomIds = category["id"].GetStrings();
            string[] atomNames = category["label_atom_id"].GetStrings();
            string[] seqIds = category["label_seq_id"].GetStrings();
            string[] compIds = category["label_comp_id"].GetStrings();
            string[] asymIds = category["label_asym_id"].GetStrings();
            string[] entityIds = category["label_entity_id"].GetStrings();           
            Console.WriteLine("row \tatom \ta.name\t  comp seq\tasym \tentity");
            foreach (int iRow in rows){
                Console.WriteLine($"{iRow}:\t {atomIds[iRow]}\t {atomNames[iRow]}\t   {compIds[iRow]}  {seqIds[iRow]} \t {asymIds[iRow]}\t {entityIds[iRow]}");
            }

            throw new NotImplementedException();
        }
    }
}
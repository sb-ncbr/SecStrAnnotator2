using System;
using System.IO;
using System.Linq;
using SecStrAnnot2.Cif;
using SecStrAnnot2.Cif.Tables;

namespace SecStrAnnot2
{
    class Program
    {
        public static DateTime SetTextDone;
        public static DateTime LexicalDone;
        public static DateTime ExtractNamesDone;
        public static DateTime SyntacticDone;
        static void Main(string[] args)
        {
            // args = new string[]{"../SecStrAnnot2_data/3ejb_updated.cif"};
            foreach (string filename in args){
                // Console.Error.WriteLine("\n" + filename);
                    try {
                        string text = "";
                        DateTime t0 = DateTime.Now;
                        using (StreamReader r = new StreamReader(filename)){
                            text = r.ReadToEnd();
                        }
                        DateTime t1 = DateTime.Now;

                        CifPackage pack = CifPackage.FromString(text);
                        string blockName = pack.BlockNames[0];
                        CifBlock block = pack[blockName];
                        
                        DateTime t2 = DateTime.Now;

                        CifCategory category = block["_atom_site"];
                        // CifItem item = category["label_atom_id"];
                        // int[] cAlphas = block.GetItem("_atom_site.label_atom_id").GetRowsWith("CA"); // or item.GetRowsWith("CA");
                        // AtomTable atoms = new AtomTable(block, cAlphas); // or new AtomTable(category, cAlphas);                        

                        // if (block.ContainsCategory("_pdbx_struct_sheet_hbond")){
                        //     CifCategory sheetCategory = block["_pdbx_struct_sheet_hbond"];
                        //     CifItem sheetId = sheetCategory["sheet_id"];
                        //     int[] sheetIdB = sheetId.GetRowsWith("B");
                        //     Table table = sheetCategory.MakeTable(sheetIdB, ("sheet_id", CifValueType.Char), ("range_1_label_comp_id", CifValueType.String), ("range_1_label_seq_id", CifValueType.Integer));
                        // }

                        // int[] startsOfEntities;
                        // int[] startsOfChains;
                        // int[] startChainsOfEntities;
                        // int[] groupedRows = Enumerable.Range(0, category.RowCount).ToArray();
                        // groupedRows = category["label_entity_id"].GetRowsGroupedByValue(groupedRows, out startsOfEntities);
                        // Lib.LogList("Grouped rows", groupedRows);
                        // Lib.LogList("Starts of entities", startsOfEntities);
                        // groupedRows = category["label_asym_id"].GetRowsGroupedByValueInEachRegion(groupedRows, startsOfEntities, out startsOfChains, out startChainsOfEntities);
                        // Lib.LogList("Starts of chains", startsOfChains);
                        // Lib.LogList("Start chains of entities", startChainsOfEntities);
                        // return;
                        ModelCollection mc = ModelCollection.FromCifBlock(block);
                        Model model = mc.GetModel(0);
                        Lib.WriteLineDebug("model " + model.ModelNumber);

                        DateTime t3 = DateTime.Now;

                        // string[] atomIds = category["id"].GetStrings();
                        // string[] atomNames = category["label_atom_id"].GetStrings();
                        // string[] seqIds = category["label_seq_id"].GetStrings();
                        // string[] compIds = category["label_comp_id"].GetStrings();
                        // string[] asymIds = category["label_asym_id"].GetStrings();
                        // string[] entityIds = category["label_entity_id"].GetStrings();
                        // Console.WriteLine("row \tatom \ta.name\t  comp seq\tasym \tentity");
                        // for (int iRow = 0; iRow < atomIds.Length; iRow++){
                        //     Console.WriteLine($"{iRow}:\t {atomIds[iRow]}\t {atomNames[iRow]}\t   {compIds[iRow]}  {seqIds[iRow]} \t {asymIds[iRow]}\t {entityIds[iRow]}");
                        // }
                        // Lib.WriteLineDebug("SORTED:");
                        // Console.WriteLine("row \tatom \ta.name\t  comp seq\tasym \tentity");
                        // foreach (int iRow in groupedRows){
                        //     Console.WriteLine($"{iRow}:\t {atomIds[iRow]}\t {atomNames[iRow]}\t   {compIds[iRow]}  {seqIds[iRow]} \t {asymIds[iRow]}\t {entityIds[iRow]}");
                        // }

                        // Lib.WriteLineDebug(table.GetColumn<char>("sheet_id").Enumerate());
                        // Lib.WriteLineDebug(table.GetColumn<string>("range_1_label_comp_id").Enumerate());
                        // Lib.WriteLineDebug(table.GetColumn<int>("range_1_label_seq_id").Enumerate());
                        
                        //int[] cAlphas = block.GetItem("_atom_site.label_atom_id").GetIndicesWhere(name => name == "CA");
                        //int[] cAlphas = block.GetItem("_atom_site.label_atom_id").GetIndicesWhere((txt,i,j) => j-i == 2 && txt[i] == 'C' && txt[i+1] == 'A');
                        //int[] cAlphas = block.GetItem("_atom_site.label_atom_id").GetIndicesWith("CA");
                        /*double[] xs =  block.GetItem("_atom_site.Cartn_x").GetDoubles(cAlphas);
                        double[] ys =  block.GetItem("_atom_site.Cartn_y").GetDoubles(cAlphas);
                        double[] zs =  block.GetItem("_atom_site.Cartn_z").GetDoubles(cAlphas);
                        string[] names =  block.GetItem("_atom_site.auth_atom_id").GetStrings(cAlphas);*/
                        /*for (int i = 0; i < 10; i++)
                        {
                            int[] cAlphas = block.GetItem("_atom_site.label_atom_id").GetIndicesWith("CA");
                            // double[] xs =  block.GetItem("_atom_site.Cartn_x").GetDoubles();
                            // double[] ys =  block.GetItem("_atom_site.Cartn_y").GetDoubles();
                            // double[] zs =  block.GetItem("_atom_site.Cartn_z").GetDoubles();
                            // int[] resis =  block.GetItem("_atom_site.auth_seq_id").GetIntegers();
                            // Lib.WriteLineDebug(xs.Length + " xs: " + xs.Enumerate());
                            // Lib.WriteLineDebug(resis.Length + " resis: " + resis.Enumerate());
                        }*/
                        
                        //Lib.WriteLineDebug(cAlphas.Length + " CAs: " + cAlphas.Enumerate());
                        //Lib.WriteLineDebug(" resis: " + atoms.labelResSeq.Enumerate());
                        //Lib.WriteLineDebug(" Xs: " + atoms.X.Enumerate());


                        Func<TimeSpan,string> Format = span => span.TotalSeconds.ToString("0.000");

                        Lib.WriteLineDebug("Read:    " + Format(t1-t0));
                        Lib.WriteLineDebug("Parse:   " + Format(t2-t1));
                        Lib.WriteLineDebug("    Set text:   " + Format(SetTextDone-t1));
                        Lib.WriteLineDebug("    Lexical:    " + Format(LexicalDone-SetTextDone));
                        Lib.WriteLineDebug("    Names:      " + Format(ExtractNamesDone-LexicalDone));
                        Lib.WriteLineDebug("    Syntactic:  " + Format(SyntacticDone-ExtractNamesDone));
                        Lib.WriteLineDebug("    ?:          " + Format(t2-SyntacticDone));
                        Lib.WriteLineDebug("Extract: " + Format(t3-t2));
                    } catch (CifException e) {
                        Lib.WriteErrorAndExit(e.Message);
                    }
            }
        }

        private static bool CmpStr(string text, int index, string sample){
            //if (index + sample.Length > text.Length) return false;
            for (int i = 0; i < sample.Length; i++)
            {
                if (text[index+i] != sample[i]) return false;
            }
            return true;
        }
        private static bool CmpStr(string text, int index, string sample, int textLength, int sampleLength){
            if (index + sampleLength > textLength) return false;
            for (int i = 0; i < sampleLength; i++)
            {
                if (text[index+i] != sample[i]) return false;
            }
            return true;
        }
    }
}

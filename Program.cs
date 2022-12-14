using System;
using System.IO;
using System.Linq;
using System.Collections.Generic;

using Cif;
using Cif.Tables;
using Cif.Filtering;
using SecStrAnnotator2.Utils;

namespace SecStrAnnotator2
{
    class Program
    {
        static int Main(string[] args)
        {
            return protein.MainClass.Main_SecStrAnnot1(args);
        }

        static int TestingMain(string[] args)
        {
            if (args.Length == 0)
            {
                args = new string[] { "../SecStrAnnot2_data/1tqn_updated.cif" };
            }
            foreach (string filename in args)
            {
                // Console.Error.WriteLine("\n" + filename);
                Cif.Components.Protein p = CifWrapperForSecStrAnnot1.ProteinFromCifFile(filename);
                // p.Save(filename + "-converted.pdb");
                Lib2.WriteLineDebug("Read and pseudoconverted " + filename);
                // continue;

                try
                {
                    string text = "";
                    DateTime t0 = DateTime.Now;
                    using (StreamReader r = new StreamReader(filename))
                    {
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
                    ModelCollection mc = ModelCollection.FromCifBlock(block/*, block["_atom_site"]["label_atom_id"].GetRowsWith("CA")*/);
                    foreach (Model model in mc.GetModels())
                    {
                        // Lib.WriteLineDebug("Model " + model.ModelNumber + " (" + model.Atoms.Count + " atoms)\n");
                        // Lib.WriteLineDebug(model.Print() + "\n");
                    }

                    DateTime t3 = DateTime.Now;

                    Filter filter1 = Filter.IntegerInRange("label_seq_id", (100, 105), (200, 202))
                                     //  & Filter.StringEquals("label_comp_id", new string[]{"TRP","CYS"})
                                     & Filter.StringEquals("label_atom_id", new string[] { "CA", "N" })
                                     & Filter.StringEquals("label_asym_id", new string[] { "A" });
                    Filter filter2 = Filter.IntegerInRange("label_seq_id", (100, 105))
                                     | Filter.StringEquals("label_comp_id", new string[] { "TRP", "CYS" })
                                     | Filter.StringEquals("label_atom_id", new string[] { "CA" });
                    Filter filter3 = !Filter.IntegerInRange("label_seq_id", (0, 390)) & !Filter.IsNull("label_seq_id");
                    Filter filter4 = Filter.TheseRows(new int[] { 0, 1, 2, 5 });
                    Filter filter5 = Filter.Where("label_seq_id", str => str.Contains('5')) & Filter.StringEquals("label_atom_id", new string[] { "CA" });

                    int[] filteredRows = filter5.GetFilteredRows(category).ToArray();
                    string[] atomIds = category["id"].GetStrings(filteredRows);
                    string[] atomNames = category["label_atom_id"].GetStrings(filteredRows);
                    string[] seqIds = category["label_seq_id"].GetStrings(filteredRows);
                    string[] compIds = category["label_comp_id"].GetStrings(filteredRows);
                    string[] asymIds = category["label_asym_id"].GetStrings(filteredRows);
                    string[] entityIds = category["label_entity_id"].GetStrings(filteredRows);
                    Console.WriteLine("row \tatom \ta.name\t  comp seq\tasym \tentity");
                    for (int iRow = 0; iRow < atomIds.Length; iRow++)
                    {
                        Console.WriteLine($"{iRow}:\t {atomIds[iRow]}\t {atomNames[iRow]}\t   {compIds[iRow]}  {seqIds[iRow]} \t {asymIds[iRow]}\t {entityIds[iRow]}");
                    }

                    Console.WriteLine(string.Join(" ", filteredRows));
                    Console.WriteLine(category.MakeCifString(filteredRows));

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


                    Func<TimeSpan, string> Format = span => span.TotalSeconds.ToString("0.000");

                    Lib2.WriteLineDebug("Read:    " + Format(t1 - t0));
                    Lib2.WriteLineDebug("Parse:   " + Format(t2 - t1));
                    Lib2.WriteLineDebug("    Set text:   " + Format(pack.Parser.TimeStamps.SetTextDone - t1));
                    Lib2.WriteLineDebug("    Lexical:    " + Format(pack.Parser.TimeStamps.LexicalAnalysisDone - pack.Parser.TimeStamps.SetTextDone));
                    Lib2.WriteLineDebug("    Names:      " + Format(pack.Parser.TimeStamps.ExtractNamesDone - pack.Parser.TimeStamps.LexicalAnalysisDone));
                    Lib2.WriteLineDebug("    Syntactic:  " + Format(pack.Parser.TimeStamps.SyntacticAnalysisDone - pack.Parser.TimeStamps.ExtractNamesDone));
                    Lib2.WriteLineDebug("    ?:          " + Format(t2 - pack.Parser.TimeStamps.SyntacticAnalysisDone));
                    Lib2.WriteLineDebug("Extract: " + Format(t3 - t2));
                }
                catch (CifException e)
                {
                    Lib2.WriteErrorAndExit(e.Message);
                }
            }
            return 0;
        }

        private static bool CmpStr(string text, int index, string sample)
        {
            for (int i = 0; i < sample.Length; i++)
            {
                if (text[index + i] != sample[i]) return false;
            }
            return true;
        }

        private static bool CmpStr(string text, int index, string sample, int textLength, int sampleLength)
        {
            if (index + sampleLength > textLength) return false;
            for (int i = 0; i < sampleLength; i++)
            {
                if (text[index + i] != sample[i]) return false;
            }
            return true;
        }
    }
}

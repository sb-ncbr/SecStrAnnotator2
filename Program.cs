using System;
using System.IO;
using System.Linq;
using SecStrAnnot2.Cif;
using SecStrAnnot2.Cif.Tables;

namespace SecStrAnnot2
{
    class Program
    {
        static void Main(string[] args)
        {
            foreach (string filename in args){
                Console.Error.WriteLine("\n" + filename);
                    try {
                        string text = "";
                        DateTime t0 = DateTime.Now;
                        using (StreamReader r = new StreamReader(filename)){
                            text = r.ReadToEnd();
                        }
                        DateTime t1 = DateTime.Now;
                        //MyCifPack pack = MyCifPack.FromFile(filename);
                        CifPackage pack = CifPackage.FromString(text);
                        string blockName = pack.BlockNames[0];
                        CifBlock block = pack[blockName];
                        CifCategory category = block["_atom_site"];
                        //Console.WriteLine(category.Name + "[" + category.ItemKeywordNames.Length + "]\n\t" + category.ItemKeywordNames.Enumerate("\n\t"));
                        DateTime t2 = DateTime.Now;
                        AtomTable atoms = new AtomTable(block);
                        DateTime t3 = DateTime.Now;

                        Lib.WriteLineDebug("Read:    " + (t1-t0));
                        Lib.WriteLineDebug("Parse:   " + (t2-t1));
                        Lib.WriteLineDebug("Extract: " + (t3-t2));
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


using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Cif.Raw;

namespace Cif
{
    public class CifPackage
    {
        internal CifParser Parser;  // should be private?
        public string[] BlockNames { get; private set; }
        public CifBlock[] Blocks { get; private set; }
        private Dictionary<string,int> blockIndex;
        
        public CifBlock this[string blockName] => GetBlock(blockName);        
        public CifBlock GetBlock(string blockName) {
            if (!blockIndex.ContainsKey(blockName)){
                throw new KeyNotFoundExceptionWithKey<string>(blockName, "This " + this.GetType() + " does not contain block with name '" + blockName + "'");
            }
            return Blocks[blockIndex[blockName]];
        }
        public bool ContainsBlock(string blockName) => blockIndex.ContainsKey(blockName);

        public static CifPackage FromString(string text){
            return new CifPackage(text);
        }
        public static CifPackage FromFile(string filename){
            string text;
            using (StreamReader r = new StreamReader(filename)){
                text = r.ReadToEnd();
            }
            return new CifPackage(text);
        }
        private CifPackage(string text){
            this.Parser = new CifParser(text);
            this.BlockNames = Parser.BlockNames;
            this.Blocks = new CifBlock[BlockNames.Length];
            this.blockIndex = new Dictionary<string, int>(BlockNames.Length);
            for (int i = 0; i < BlockNames.Length; i++)
            {
                Blocks[i] = new CifBlock(Parser, i);
                blockIndex[BlockNames[i]] = i;
            }
        }

        

    }
}
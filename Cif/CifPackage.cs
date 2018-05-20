
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using SecStrAnnot2.Cif.Raw;

namespace SecStrAnnot2.Cif
{
    public class CifPackage
    {
        private CifParser parser;
        public string[] BlockNames { get; private set; }
        public CifBlock[] Blocks { get; private set; }
        private Dictionary<string,int> blockIndex;
        
        public CifBlock this[string blockName] => GetBlock(blockName);        
        public CifBlock GetBlock(string blockName) {
            if (!blockIndex.ContainsKey(blockName)){
                throw new KeyNotFoundException("This " + this.GetType() + " does not contain block with name '" + blockName + "'");
            }
            return Blocks[blockIndex[blockName]];
        }

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
            this.parser = new CifParser(text);
            this.BlockNames = parser.BlockNames;
            this.Blocks = new CifBlock[BlockNames.Length];
            this.blockIndex = new Dictionary<string, int>(BlockNames.Length);
            for (int i = 0; i < BlockNames.Length; i++)
            {
                Blocks[i] = new CifBlock(parser, i);
                blockIndex[BlockNames[i]] = i;
            }
        }

        

    }
}
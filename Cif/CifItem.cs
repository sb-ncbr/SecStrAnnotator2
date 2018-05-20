
using System;
using System.Collections.Generic;
using System.Linq;
using SecStrAnnot2.Cif.Raw;

namespace SecStrAnnot2.Cif
{
    public class MyCifItem
    {
        private CifParser parser;
        private int iTag;
        public string KeywordName { get; private set; }
        public string FullName { get; private set; }
        public int Count { get; private set; }

        internal MyCifItem(CifParser parser, string keyword, int iTag){
            this.parser = parser;
            this.iTag = iTag;
            this.FullName = parser.TagNames[iTag];
            this.KeywordName = FullName.Split('.', 2)[1];
            this.Count = parser.CountValuesForTag(iTag);
        }
        
        public object GetValues(CifValueType type){
            switch(type){
                case CifValueType.String:
                    return GetStrings();
                case CifValueType.Char:
                    return GetChars();
                case CifValueType.Integer:
                    return GetIntegers();
                case CifValueType.Double:
                    return GetDoubles();
                default:
                    throw new NotImplementedException();
            }
        }

        public string[] GetStrings (){
            return parser.GetValuesAsStrings(iTag);
        }
        public char[] GetChars (){
            return parser.GetValuesAsChars(iTag);
        }
        public int[] GetIntegers (){
            return parser.GetValuesAsIntegers(iTag);
        }
        public double[] GetDoubles (){
            return parser.GetValuesAsDoubles(iTag);
        }

        //TODO Add methods for obtaining only a subset of values.
    }
}
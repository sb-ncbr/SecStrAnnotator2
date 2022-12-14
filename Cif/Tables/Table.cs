using System;
using System.Collections.Generic;
using System.Linq;
using Cif.Libraries;

namespace Cif.Tables
{
    public class Table
    {
        public int RowCount { get; protected set; }
        public int ColumnCount { get; protected set; }
        public string[] ColumnNames { get; protected set; }
        protected object[] columnArray;
        private CifValueType[] cifTypeArray;
        private Type[] backtypeArray;
        protected Dictionary<string, int> columnIndex;

        private static Dictionary<CifValueType,Type> typeDictionary = new Dictionary<CifValueType,Type>{ 
            { CifValueType.String, "".GetType() },
            { CifValueType.Char, 'a'.GetType() },
            { CifValueType.Integer, ((int)0).GetType() },
            { CifValueType.Double, ((double)0).GetType() }
        };

        public bool ContainsColumn(string columnName){
            return columnIndex.ContainsKey(columnName);
        }

        public CifValueType GetCifValueType(string keyword){
            if (!ContainsColumn(keyword)){
                throw new KeyNotFoundExceptionWithKey<string>(keyword);
            }
            return cifTypeArray[columnIndex[keyword]];
        }

        public Type GetBackingType(string keyword){
            if (!ContainsColumn(keyword)){
                throw new KeyNotFoundExceptionWithKey<string>(keyword);
            }
            return backtypeArray[columnIndex[keyword]];
        }

        public T[] GetColumn<T>(string keyword){
            if (!ContainsColumn(keyword)){
                throw new KeyNotFoundExceptionWithKey<string>(keyword);
            }
            return columnArray[columnIndex[keyword]] as T[];
        }


        public Table (CifCategory cifCategory, params (string, CifValueType)[] schema) { 
            ColumnCount = schema.Length;
            ColumnNames = new string[ColumnCount];
            columnArray = new object[ColumnCount];
            cifTypeArray = new CifValueType[ColumnCount];
            backtypeArray = new Type[ColumnCount];
            columnIndex = new Dictionary<string, int>(ColumnCount);

            for (int i = 0; i < ColumnCount; i++)
            {
                (string keyword, CifValueType cifType) = schema[i];
                ColumnNames[i] = keyword;
                if (cifCategory.ContainsItem(keyword)){
                    columnArray[i] = cifCategory[keyword].GetValues(cifType);
                } else {
                    throw new KeyNotFoundExceptionWithKey<string>(keyword, "CIF category " + cifCategory.Name + " does not contain item with keyword " + keyword);
                }
                cifTypeArray[i] = cifType;
                backtypeArray[i] = typeDictionary[cifType];
                columnIndex[keyword] = i; 
            }
        }

        public Table (CifCategory cifCategory, int[] rows, params (string, CifValueType)[] schema) { 
            ColumnCount = schema.Length;
            ColumnNames = new string[ColumnCount];
            columnArray = new object[ColumnCount];
            cifTypeArray = new CifValueType[ColumnCount];
            backtypeArray = new Type[ColumnCount];
            columnIndex = new Dictionary<string, int>(ColumnCount);

            for (int i = 0; i < ColumnCount; i++)
            {
                (string keyword, CifValueType cifType) = schema[i];
                ColumnNames[i] = keyword;
                if (cifCategory.ContainsItem(keyword)){
                    columnArray[i] = cifCategory[keyword].GetValues(cifType, rows);
                } else {
                    throw new KeyNotFoundExceptionWithKey<string>(keyword, "CIF category " + cifCategory.Name + " does not contain item with keyword " + keyword);
                }
                cifTypeArray[i] = cifType;
                backtypeArray[i] = typeDictionary[cifType];
                columnIndex[keyword] = i; 
            }
        }

        public void Print(){
            for (int i = 0; i < ColumnCount; i++)
            {
                Console.WriteLine(ColumnNames[i] + ":");
                Console.WriteLine(string.Join (", ", columnArray[i] as IEnumerable<object>));
            }
        }
    }
}
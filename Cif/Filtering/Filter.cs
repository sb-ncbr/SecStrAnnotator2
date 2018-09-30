using System;
using System.Collections.Generic;
using System.Linq;
using Cif.Raw;

namespace Cif.Filtering
{
    public abstract class Filter
    {
        public virtual List<int> GetFilteredRows(CifCategory category){
            List<int> allRows = Enumerable.Range(0, category.RowCount).ToList();  // Could be more efficient?
            return GetFilteredRows(category, allRows);
        }
        protected abstract List<int> GetFilteredRows(CifCategory category, List<int> inputRows);


        private abstract class OneColumnFilter : Filter {
            private string columnName;
            protected OneColumnFilter(string columnName){
                this.columnName = columnName;
            }
            protected (CifParser, int) GetParserAndITag(CifCategory category){
                return (category.parser, category.GetItem(columnName).iTag);
            }
        }


        private class StringEqualsFilter : OneColumnFilter {
            private string[] values;
            public StringEqualsFilter(string columnName, string[] values) : base(columnName) {
                this.values = values;
            }
            protected override List<int> GetFilteredRows(CifCategory category, List<int> inputRows){
                (CifParser parser, int iTag) = GetParserAndITag(category);
                if (values.Length == 0) {
                    return new List<int>();
                } else if (values.Length == 1) {
                    string value  = values[0];
                    return parser.GetIndicesWith(iTag, value, inputRows);
                } else {
                    return parser.GetIndicesWith(iTag, values, inputRows);
                }
            }
        }


        private class PredicateFilter : OneColumnFilter {
            private Func<string,bool> predicate;
            public PredicateFilter(string columnName, Func<string,bool> predicate) : base(columnName) {
                this.predicate = predicate;
            }
            protected override List<int> GetFilteredRows(CifCategory category, List<int> inputRows){
                (CifParser parser, int iTag) = GetParserAndITag(category);
                return parser.GetIndicesWhere(iTag, predicate, inputRows);
            }
        }


        private class IntegerInRangeFilter : OneColumnFilter {
            private (int, int)[] ranges;
            public IntegerInRangeFilter(string columnName, (int,int)[] ranges) : base(columnName) {
                this.ranges = ranges;
            }
            protected override List<int> GetFilteredRows(CifCategory category, List<int> inputRows){
                (CifParser parser, int iTag) = GetParserAndITag(category);
                if (ranges.Length == 0) {
                    return new List<int>();
                } else if (ranges.Length == 1) {
                    (int,int) range  = ranges[0];
                    return parser.GetIndicesWithIntegerInRange(iTag, range, inputRows);
                } else {
                    return parser.GetIndicesWithIntegerInRange(iTag, ranges, inputRows);
                }
            }
        }


        private class TheseRowsFilter : Filter {
            private List<int> rows;
            public TheseRowsFilter(List<int> rows) {
                this.rows = rows;
            }
            protected override List<int> GetFilteredRows(CifCategory category, List<int> inputRows){
                return Intersect(inputRows, this.rows);
            }
            private List<int> Intersect(List<int> first, List<int> second){
                int n1 = first.Count;
                int n2 = second.Count;
                List<int> result = new List<int>(Math.Min(n1, n2));
                if (n1 == 0 || n2 == 0){
                    return result;
                } else {
                    int i = 0;
                    int j = 0;
                    while (i < n1 && j < n2) {
                        int val1 = first[i];
                        int val2 = second[j];
                        if (val1 < val2) {
                            // don't add nothing
                            i++;
                        } else if (val2 < val1) {
                            // don't add nothing
                            j++;
                        } else { // val1 == val2
                            result.Add(val1);
                            i++;
                            j++;
                        }
                    }
                    return result;
                }
            }
        }


        private class AndFilter : Filter {
            private Filter[] filters;
            public AndFilter(Filter[] filters){
                this.filters = filters;
            }
            protected override List<int> GetFilteredRows(CifCategory category, List<int> inputRows){
                List<int> rows = inputRows;
                foreach(Filter filter in this.filters){
                    rows = filter.GetFilteredRows(category, rows);
                }
                return rows;
            }
        }


        private class OrFilter : Filter {
            private Filter[] filters;
            public OrFilter(Filter[] filters){
                this.filters = filters;
            }
            protected override List<int> GetFilteredRows(CifCategory category, List<int> inputRows){
                List<int> outputRows = new List<int>();
                foreach(Filter filter in this.filters){
                    List<int> newRows = filter.GetFilteredRows(category, inputRows);
                    outputRows = Merge(outputRows, newRows);
                }
                return outputRows;
            }
            private List<int> Merge(List<int> first, List<int> second){
                int n1 = first.Count;
                int n2 = second.Count;
                if (n1 == 0){
                    return second;
                } else if (n2 == 0){
                    return first;
                } else {
                    List<int> result = new List<int>(Math.Max(n1, n2));
                    int i = 0;
                    int j = 0;
                    while (i < n1 && j < n2) {
                        int val1 = first[i];
                        int val2 = second[j];
                        if (val1 < val2) {
                            result.Add(val1);
                            i++;
                        } else if (val2 < val1) {
                            result.Add(val2);
                            j++;
                        } else { // val1 == val2
                            result.Add(val1);
                            i++;
                            j++;
                        }
                    }
                    if (i < n1) {
                        for (; i < n1; i++)
                        result.Add(first[i]);
                    } else if (j < n2) {
                        for (; j < n2; j++)
                        result.Add(second[j]);
                    }
                    return result;
                }
            }
        }


        private class NotFilter : Filter {
            private Filter filter;
            public NotFilter(Filter filter){
                this.filter = filter;
            }
            protected override List<int> GetFilteredRows(CifCategory category, List<int> inputRows){
                List<int> omitted = filter.GetFilteredRows(category, inputRows);
                return SetMinus(inputRows, omitted);
            }
            private List<int> SetMinus(List<int> first, List<int> second){
                int n1 = first.Count;
                int n2 = second.Count;
                if (n1 == 0){
                    return new List<int>();
                } else if (n2 == 0){
                    return first;
                } else {
                    List<int> result = new List<int>(Math.Max(n1, n2));
                    int i = 0;
                    int j = 0;
                    while (i < n1 && j < n2) {
                        int val1 = first[i];
                        int val2 = second[j];
                        if (val1 < val2) {
                            result.Add(val1);
                            i++;
                        } else if (val1 == val2) {
                            // do not add
                            i++;
                            j++;
                        } else { // val2 < val1
                            // do nothing (missing value in first set)
                            j++;
                        }
                    }
                    if (i < n1) {
                        for (; i < n1; i++)
                        result.Add(first[i]);
                    }
                    return result;
                }
            }
        
        }


        public static Filter IsNull(string columnName) => new StringEqualsFilter(columnName, new string[] {".", "?"});
        public static Filter IsDot(string columnName) => new StringEqualsFilter(columnName, new string[] {"."});
        public static Filter IsQuestionMark(string columnName) => new StringEqualsFilter(columnName, new string[] {"?"});

        public static Filter StringEquals(string columnName, string value) => new StringEqualsFilter(columnName, new string[] {value});
        public static Filter StringEquals(string columnName, string[] oneOfValues) => new StringEqualsFilter(columnName, oneOfValues);

        public static Filter IntegerInRange(string columnName, (int,int) range) => new IntegerInRangeFilter(columnName, new (int,int)[] {range});
        public static Filter IntegerInRange(string columnName, params (int,int)[] ranges) => new IntegerInRangeFilter(columnName, ranges);

        public static Filter Where(string columnName, Func<string,bool> predicate) => new PredicateFilter(columnName, predicate);
        
        public static Filter TheseRows(IEnumerable<int> rows) => new TheseRowsFilter(rows.ToList());

        public static Filter And(params Filter[] filters) => new AndFilter(filters);
        public static Filter Or(params Filter[] filters) => new OrFilter(filters);
        public static Filter Not(Filter filter) => new NotFilter(filter);

        public static Filter operator & (Filter first, Filter second) => And(first, second);
        public static Filter operator | (Filter first, Filter second) => Or(first, second);
        public static Filter operator ! (Filter filter) => Not(filter);
    }
}
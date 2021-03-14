using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;
using System.Reflection;

namespace protein.Libraries {
    static class Lib {
        private const ConsoleColor ERROR_COLOR = ConsoleColor.Red;
        private const ConsoleColor WARNING_COLOR = ConsoleColor.Red; // ConsoleColor.Magenta;
        private const ConsoleColor DEBUG_COLOR = ConsoleColor.Blue;
        public static bool DoWriteDebug = false;

        /** Reads a CSV-formatted table from StreamReader r. 
			The table must have numColumns columns, lines are separated by '\n', fields in lines are separated by separator.
			If commentSign is not null, everything that is found on a line after commentSign is ignored. 
			Returns a list of lines, each line is a list of fields. */
        public static List<List<String>> ReadCSV(StreamReader r, int numColumns, char separator, char? commentSign) {
            List<List<String>> result = new List<List<String>>();
            String line;
            while ((line = r.ReadLine()) != null) {
                String[] fields = commentSign != null ? line.Split(new char[] { (char)commentSign }, 2)[0].Split(separator) : line.Split(separator);
                int lenField = fields.Length;
                while (lenField > 0 && (fields[lenField - 1].Trim() == ""))
                    lenField--;
                if (lenField == numColumns)
                    result.Add(fields.ToList().GetRange(0, numColumns));
                else if (lenField != 0)
                    throw new FormatException("ReadCSV: Error on line \"" + line + ".");
            }
            return result;
        }

        /** Writes a CSV-formatted table to StreamWriter w. 
			Each line of the table must have the same number of fields, i.e. each item in the List table is a List with the same number of items.
			Lines in output are separated by '\n', fields in lines are separated by separator. */
        public static void WriteCSV(StreamWriter w, List<List<String>> table, char separator) {
            HashSet<int> lengths = new HashSet<int>(table.Select(x => x.Count));
            if (lengths.Count == 0) {
                // empty table, no need to print anything :)
            } else if (lengths.Count == 1) {
                // all lines have the same number of fields, OK
                foreach (List<String> line in table) {
                    for (int i = 0; i <= line.Count - 2; i++)
                        w.Write(line[i] + separator);
                    if (line.Count > 0)
                        w.WriteLine(line.Last());
                }
            } else if (lengths.Count > 1) {
                // not all lines have the same number of fields, not OK
                throw new InvalidDataException("WriteCSV: not all lines in the table have the same number of fields.");
            }
        }

        /** Does the same as WriteCSV() but treats the last column of the table as the comments, i.e. prepends the commentSign to them.
		 * The parameter introComment (may be multi-line) will be printed as a comment in the beginning of the file.
			Each line of the table must have the same number of fields, i.e. each item in the List table is a List with the same number of items.
			Lines in output are separated by '\n', fields in lines are separated by separator. */
        public static void WriteCSVWithComments(StreamWriter w, List<List<String>> table, char separator, String commentSign, String introComment) {
            if (introComment != null)
                w.WriteLine(String.Concat(introComment.Split('\n').Select(x => x != "" ? (commentSign + x + "\n") : "\n")));
            Lib.WriteCSV(w, table.Select(line => line.Select((x, i) => i == line.Count - 1 && x != "" ? (commentSign + x) : x).ToList()).ToList(), separator);
        }


        /** Runs a command with specified arguments. */
        public static bool RunCommand(String command, String args) {
            bool rooted = Path.IsPathRooted(command);
            if (rooted) {
                if (!File.Exists(command)) {
                    Lib.WriteError("File '{0}' not found.", command);
                    return false;
                }
            } else {
                String commandHere = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, command);
                if (File.Exists(commandHere)) {
                    Console.WriteLine($"Info: File '{command}; found in the base directory '{AppDomain.CurrentDomain.BaseDirectory}'");
                    command = commandHere;
                } else {
                    Console.WriteLine($"Info: File '{command}' not found in the base directory '{AppDomain.CurrentDomain.BaseDirectory}'.");
                    Console.WriteLine($"      Running '{command}' from current directory or $PATH.");
                }
            }
            try {
                Process proc = new Process();
                proc.StartInfo = new ProcessStartInfo();
                // proc.StartInfo.RedirectStandardInput = false;
                // proc.StartInfo.RedirectStandardOutput = false;
                // proc.StartInfo.RedirectStandardError = false;
                proc.StartInfo.FileName = command;
                proc.StartInfo.Arguments = args;
                proc.StartInfo.UseShellExecute = false; // without this it throws something about xdg-open on Linux if the pymol is not installed
                Console.WriteLine("Executing: ");
                Console.WriteLine("    " + command + " " + args);
                proc.Start();
                proc.WaitForExit();
                int exitCode = proc.ExitCode;
                bool success = exitCode == 0;
                String successString = success ? "success" : "fail";
                Console.WriteLine($"Exit code: {exitCode} ({successString})\n");
                return success;
            } catch (System.ComponentModel.Win32Exception e) {
                Lib.WriteError($"Could not run: {command}");
                Lib.WriteError($"Win32Exception: {e.Message}\n");
                return false;
            }
        }


        /** Runs PyMOL giving it a python script from file scriptToBeRun. 
		* This means running: pymol -qcyr scriptToBeRun.*/
        public static bool RunPyMOLScriptWithCommandLineArguments(String pymolExecutable, String scriptToBeRun, IEnumerable<String> arguments) {
            if (!Path.IsPathRooted(scriptToBeRun)) {
                scriptToBeRun = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, scriptToBeRun);
            }
            if (!File.Exists(scriptToBeRun)) {
                Lib.WriteErrorAndExit("Script file '{0}' not found.", scriptToBeRun);
            }
            bool success = RunCommand(pymolExecutable, "-qcyr " + scriptToBeRun + " -- " + arguments.EnumerateWithSeparators(" "));
            return success;
        }

        /** Runs DSSP program. */
        public static bool RunDSSP(String dsspExecutable, String input, String output) {
            bool success = RunCommand(dsspExecutable, "-i " + input + " -o " + output);
            return success;
        }

        public static void WriteInColor(ConsoleColor color, String s, params object[] args) {
            ConsoleColor orig = Console.ForegroundColor;
            Console.ForegroundColor = color;
            Console.Write(String.Format(s, args));
            Console.ForegroundColor = orig;
        }

        public static void WriteDebug(String s, params object[] args) {
            if (DoWriteDebug) {
                ConsoleColor orig = Console.ForegroundColor;
                Console.ForegroundColor = DEBUG_COLOR;
                Console.Write(String.Format(s, args));
                Console.ForegroundColor = orig;
            }
        }

        public static void WriteLineDebug(String s, params object[] args) {
            if (DoWriteDebug) {
                ConsoleColor orig = Console.ForegroundColor;
                Console.ForegroundColor = DEBUG_COLOR;
                Console.WriteLine(String.Format(s, args));
                Console.ForegroundColor = orig;
            }
        }
        public static void WriteLineDebug(string s) => WriteLineDebug(s, new object[] { });

        public static void WriteError(String s, params object[] args) {
            ConsoleColor orig = Console.ForegroundColor;
            Console.ForegroundColor = ERROR_COLOR;
            Console.Error.WriteLine("ERROR: " + String.Format(s, args));
            Console.ForegroundColor = orig;
        }

        public static void WriteErrorAndExit(String s, params object[] args) {
            WriteError(s, args);
            Environment.Exit(-1);
        }

        public static void WriteWarning(String s, params object[] args) {
            ConsoleColor orig = Console.ForegroundColor;
            Console.ForegroundColor = WARNING_COLOR;
            Console.Error.WriteLine("WARNING: " + String.Format(s, args));
            Console.ForegroundColor = orig;
        }

        public static Version BuildVersion {
            get {
                return TimeToBuildVersion(RealBuildTime);
                return Assembly.GetExecutingAssembly().GetName().Version;
            }
        }

        public static DateTime RealBuildTime {
            get {
                var assembly = typeof(Lib).GetTypeInfo().Assembly;
                string assemblyName = assembly.GetName().Name;
                Stream resource = assembly.GetManifestResourceStream($"{assemblyName}.buildtime");
                string time;
                using (StreamReader r = new StreamReader(resource)) {
                    time = r.ReadToEnd();
                }
                resource.Close();
                return DateTime.Parse(time);
            }
        }

        public static DateTime BuildTime {
            get {
                return RealBuildTime;
                // Version version = BuildVersion;
                // TimeSpan span = new TimeSpan (
                // 	               TimeSpan.TicksPerDay * version.Build + // days since 1 January 2000
                // 	               TimeSpan.TicksPerSecond * 2 * version.Revision);
                // DateTime buildTime = new DateTime (2000, 1, 1).Add (span);
                // return buildTime;
            }
        }

        private static Version TimeToBuildVersion(DateTime time) {
            TimeSpan since2000 = time.Date - new DateTime(2000, 1, 1);
            return new Version(0, 0, (int)since2000.TotalDays, (int)time.TimeOfDay.TotalSeconds % 100000);
        }

        /**Orders the IEnumerable in the same way as OrderBy() and returns the ranks of the elements of the original IEnumerable.
		 * e.g. {a,d,b,c} --> result={a,b,c,d}, ranks={1,4,2,3}*/
        public static IEnumerable<TSource> OrderAndGetRanks<TSource, TKey>(IEnumerable<TSource> list, Func<TSource, TKey> keySelector, out int[] ranks) {
            IEnumerable<(TSource, int)> ordered =
                list.Select((x, origIndex) => (x, origIndex))
                    .OrderBy(t => keySelector(t.Item1));
            ranks =
                ordered
                    .Select((t, orderedIndex) => (t.Item2, orderedIndex))
                    .OrderBy(t => t.Item1)
                    .Select(t => t.Item2)
                    .ToArray();
            return
                ordered.Select(t => t.Item1);
        }

        public class Shuffler {
            private int[] oldToNewIndex;
            private int[] newToOldIndex;
            public int MaxNewIndex => newToOldIndex.Length - 1;
            public int MaxOldIndex => oldToNewIndex.Length - 1;

            private Shuffler() { }
            public Shuffler(IEnumerable<int> oldIndices) {
                if (oldIndices.Count() > 0) {
                    this.newToOldIndex = oldIndices.ToArray();
                    int maxOld = this.newToOldIndex.Max();
                    this.oldToNewIndex = Enumerable.Range(0, maxOld + 1).Select(x => -1).ToArray();
                    for (int i = 0; i < this.newToOldIndex.Length; i++) {
                        this.oldToNewIndex[this.newToOldIndex[i]] = i;
                    }
                } else {
                    this.oldToNewIndex = new int[0];
                    this.newToOldIndex = new int[0];
                }
            }

            public Shuffler Inverted() {
                Shuffler m = new Shuffler();
                m.oldToNewIndex = this.newToOldIndex;
                m.newToOldIndex = this.oldToNewIndex;
                return m;
            }

            public static Shuffler FromMatching(IEnumerable<(int, int)> matching) {
                if (matching.Count() > 0) {
                    int maxOld = matching.Select(t => t.Item1).Max();
                    int maxNew = matching.Select(t => t.Item2).Max();
                    Shuffler shuffler = new Shuffler();
                    shuffler.oldToNewIndex = Enumerable.Repeat(-1, maxOld + 1).ToArray();
                    shuffler.newToOldIndex = Enumerable.Repeat(-1, maxNew + 1).ToArray();
                    foreach (var match in matching) {
                        shuffler.oldToNewIndex[match.Item1] = match.Item2;
                        shuffler.newToOldIndex[match.Item2] = match.Item1;
                    }
                    return shuffler;
                } else {
                    return new Shuffler(new int[] { });
                }
            }

            public bool HasNewIndex(int i) {
                return i >= 0 && i < oldToNewIndex.Length && oldToNewIndex[i] >= 0;
            }
            public bool HasOldIndex(int j) {
                return j >= 0 && j < newToOldIndex.Length && newToOldIndex[j] >= 0;
            }

            public int NewIndex(int i) {
                if (!HasNewIndex(i))
                    throw new Exception("Invalid index " + i + ".");
                return oldToNewIndex[i];
            }

            public int OldIndex(int j) {
                if (!HasOldIndex(j))
                    throw new Exception("Invalid index " + j + ".");
                return newToOldIndex[j];
            }

            public IEnumerable<T> Shuffle<T>(IEnumerable<T> elements) {
                if (newToOldIndex.Any(i => i < 0))
                    throw new Exception("Some elements of the shuffled array are missing.");
                T[] elemArray = elements.ToArray();
                return newToOldIndex.Select(j => elemArray[j]);
            }

            public IEnumerable<T> ShuffleBack<T>(IEnumerable<T> elements) {
                if (oldToNewIndex.Any(j => j < 0))
                    throw new Exception("Some elements of the original array are missing.");
                T[] elemArray = elements.ToArray();
                return oldToNewIndex.Select(j => elemArray[j]);
            }

            public IEnumerable<T> ShuffleBack<T>(IEnumerable<T> elements, Func<T> createDefaultElement) {
                T[] elemArray = elements.ToArray();
                return oldToNewIndex.Select(j => j >= 0 ? elemArray[j] : createDefaultElement());
            }

            public T[,] ShuffleRows<T>(T[,] array) {
                if (array == null)
                    return null;
                if (newToOldIndex.Any(i => i < 0))
                    throw new Exception("Some elements of the shuffled array are missing.");
                int rows = MaxNewIndex + 1;
                int columns = array.GetLength(1);
                T[,] result = new T[rows, columns];
                for (int i = 0; i < rows; i++) {
                    for (int j = 0; j < columns; j++) {
                        result[i, j] = array[newToOldIndex[i], j];
                    }
                }
                return result;
            }
            public T[,] ShuffleColumns<T>(T[,] array) {
                if (array == null)
                    return null;
                if (newToOldIndex.Any(i => i < 0))
                    throw new Exception("Some elements of the shuffled array are missing.");
                int rows = array.GetLength(0);
                int columns = MaxNewIndex + 1;
                T[,] result = new T[rows, columns];
                for (int i = 0; i < rows; i++) {
                    for (int j = 0; j < columns; j++) {
                        result[i, j] = array[i, newToOldIndex[j]];
                    }
                }
                return result;
            }
            public T[,] ShuffleRowsAndColumns<T>(T[,] array) {
                return ShuffleRows(ShuffleColumns(array));
            }

            public IEnumerable<(int, int)> UpdateIndices(IEnumerable<(int, int)> matching) {
                return matching.Where(t => HasNewIndex(t.Item1) && HasNewIndex(t.Item2)).Select(t => (NewIndex(t.Item1), NewIndex(t.Item2)));
            }

            public IEnumerable<(int, int, T)> UpdateIndices<T>(IEnumerable<(int, int, T)> matchingWithLabels) {
                return matchingWithLabels.Where(t => HasNewIndex(t.Item1) && HasNewIndex(t.Item2)).Select(t => (NewIndex(t.Item1), NewIndex(t.Item2), t.Item3));
            }

            public override string ToString() {
                return "Shuffler{" + oldToNewIndex.Select((x, i) => x >= 0 ? i + "->" + x : null).Where(s => s != null).EnumerateWithCommas() + "}";
            }
        }

        public static IEnumerable<TSource> OrderAndGetShuffler<TSource, TKey>(this IEnumerable<TSource> list, Func<TSource, TKey> keySelector, out Shuffler shuffler) {
            IEnumerable<(TSource, int)> ordered =
                list.Select((x, origIndex) => (x, origIndex))
                    .OrderBy(t => keySelector(t.Item1));
            shuffler = new Shuffler(ordered.Select(t => t.Item2));
            return
                ordered.Select(t => t.Item1);
        }

        public static IEnumerable<TSource> OrderAndGetShuffler<TSource>(this IEnumerable<TSource> list, out Shuffler shuffler) {
            return list.OrderAndGetShuffler(x => x, out shuffler);
        }

        public static IEnumerable<TSource> ConcatAndGetShuffler<TSource>(this IEnumerable<TSource> list, IEnumerable<TSource> second, out Shuffler shuffler) {
            int n1 = list.Count();
            int n2 = second.Count();
            shuffler = new Shuffler(Enumerable.Range(0, n2).Select(i => n1 + i)).Inverted();
            return list.Concat(second);
        }

        public static IEnumerable<T> WhereAndGetShuffler<T>(this IEnumerable<T> list, Func<T, bool> predicate, out Shuffler shuffler) {
            IEnumerable<int> indices = list.IndicesWhere(predicate);
            shuffler = new Shuffler(indices);
            return shuffler.Shuffle(list);
        }

        public static List<T> Merge<T>(List<List<T>> lists, IComparer<T> comparer) {
            List<T> result = new List<T>();
            List<T> next;
            do {
                next = null;
                foreach (List<T> list in lists) {
                    if (list.Count != 0) {
                        if (next == null || comparer.Compare(list.First(), next.First()) < 0) {
                            next = list;
                        }
                    }
                }
                if (next != null) {
                    result.Add(next.First());
                    next.RemoveAt(0);
                }
            } while (next != null);
            return result;
        }

        public static void DoForAllFiles(String directory, Action<String> action) {
            if (!Directory.Exists(directory)) return;
            foreach (String file in Directory.GetFiles(directory)) {
                action(file);
            }
            foreach (String dir in Directory.GetDirectories(directory)) {
                DoForAllFiles(dir, action);
            }
        }

        public static void DoForAllFiles(String directory, Action<String> action, String regex) {
            DoForAllFiles(directory, delegate (String file) {
                String fileName = Path.GetFileName(file);
                if (Regex.IsMatch(fileName, regex)) action(file);
            });
        }

        public static String EnumerateWithSeparators<T>(this IEnumerable<T> list, String separator) {
            return string.Join(separator, list);
            // if (list.Count() == 0)
            //     return "";
            // else
            //     return list.Skip(1).Aggregate(list.First().ToString(), (s, x) => s + separator + x.ToString());
        }

        public static String EnumerateWithCommas<T>(this IEnumerable<T> list) {
            return string.Join(", ", list);
            // return EnumerateWithSeparators(list, ", ");
        }

        public static IEnumerable<int> IndicesWhere<T>(this IEnumerable<T> list, Func<T, bool> predicate) {
            return list.Select((x, i) => predicate(x) ? i : -1).Where(i => i >= 0);
        }

        public static int ArgMin<T>(this IEnumerable<T> source, out T minValue) where T : IComparable<T> {
            if (source == null)
                throw new ArgumentNullException();
            if (source.Count() == 0)
                throw new InvalidOperationException();
            int argmin = 0;
            minValue = source.First();
            int i = 1;
            foreach (T x in source.Skip(1)) {
                if (x.CompareTo(minValue) < 0) {
                    minValue = x;
                    argmin = i;
                }
                i++;
            }
            return argmin;
        }

        public static int ArgMin<T>(this IEnumerable<T> source) where T : IComparable<T> {
            T dump;
            return source.ArgMin(out dump);
        }

        public static int ArgMax<T>(this IEnumerable<T> source, out T maxValue) where T : IComparable<T> {
            if (source == null)
                throw new ArgumentNullException();
            if (source.Count() == 0)
                throw new InvalidOperationException();
            int argmax = 0;
            maxValue = source.First();
            int i = 1;
            foreach (T x in source.Skip(1)) {
                if (x.CompareTo(maxValue) > 0) {
                    maxValue = x;
                    argmax = i;
                }
                i++;
            }
            return argmax;
        }

        public static int ArgMax<T>(this IEnumerable<T> source) where T : IComparable<T> {
            T dump;
            return source.ArgMax(out dump);
        }

        public static int? GetOrAssignNext<K>(this IDictionary<K, int?> dictionary, K key, ref int valueCounter) {
            if (dictionary.ContainsKey(key))
                return dictionary[key];
            else {
                dictionary[key] = valueCounter;
                return valueCounter++;
            }
        }

        public static void MultidictionaryAdd<K, V>(this IDictionary<K, List<V>> dictionary, K key, V value) {
            if (!dictionary.ContainsKey(key))
                dictionary[key] = new List<V>();
            dictionary[key].Add(value);
        }

        public static bool IsSorted<T>(this IEnumerable<T> list, Comparison<T> comparison) {
            if (list.Count() == 0)
                return true;
            T previous = list.First();
            foreach (T next in list.Skip(1)) {
                if (comparison(previous, next) > 0)
                    return false;
                previous = next;
            }
            return true;
        }

        public static bool IsSorted<T>(this IEnumerable<T> list) where T : IComparable<T> {
            return list.IsSorted((x, y) => x.CompareTo(y));
        }

        /**Creates a new 1D array that contains the same elements on the same positions as in original.*/
        public static T[] Resized<T>(this T[] original, int length) {
            T[] resized = new T[length];
            Array.Copy(original, resized, Math.Min(original.Length, length));
            return resized;
        }

        /**Creates a new 2D array that contains the same elements on the same positions as in original.*/
        public static T[,] Resized2D<T>(this T[,] original, int rows, int cols) {
            T[,] resized = new T[rows, cols];
            int minX = Math.Min(original.GetLength(0), resized.GetLength(0));
            int minY = Math.Min(original.GetLength(1), resized.GetLength(1));

            for (int i = 0; i < minX; ++i)
                Array.Copy(original, i * original.GetLength(1), resized, i * resized.GetLength(1), minY);

            return resized;
        }

        /**Creates a new 2D array that contains the same elements on the same positions as in original.*/
        public static T[,] Transposed<T>(this T[,] original) {
            int rows = original.GetLength(1);
            int cols = original.GetLength(0);
            T[,] transposed = new T[rows, cols];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++) {
                    transposed[i, j] = original[j, i];
                }
            return transposed;
        }

        /**Fills 2D array with values (array[i,j] = getInitValue(i,j)) and returns the array itself.*/
        public static T[,] Fill<T>(this T[,] array, Func<int, int, T> getInitValue) {
            for (int i = 0; i < array.GetLength(0); i++) {
                for (int j = 0; j < array.GetLength(1); j++) {
                    array[i, j] = getInitValue(i, j);
                }
            }
            return array;
        }

        /**Fills 2D array with values (array[i,j] = initValue) and returns the array itself.*/
        public static T[,] Fill<T>(this T[,] array, T initValue) {
            int m = array.GetLength(0);
            int n = array.GetLength(1);
            if (m > 0) {
                for (int j = 0; j < n; j++) {
                    array[0, j] = initValue;
                }
                for (int i = 1; i < m; i++) {
                    Array.Copy(array, 0, array, i * n, n);
                }
            }
            return array;
        }

        public static T[,] Select2D<S, T>(this S[,] array, Func<S, T> selector) {
            int m = array.GetLength(0);
            int n = array.GetLength(1);
            T[,] result = new T[m, n];
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                    result[i, j] = selector(array[i, j]);
                }
            }
            return result;
        }

        public static IEnumerable<T> GetRange<T>(this T[] array, int start, int count) {
            return Enumerable.Range(start, count).Select(i => array[i]);
        }

        public static bool InRanges(this IEnumerable<(int, int)> ranges, int value) {
            foreach (var range in ranges) {
                if (range.Item1 <= value && value <= range.Item2)
                    return true;
            }
            return false;
        }

        public static string Repeat(this string str, int repeats) {
            if (str.Length == 0 || repeats == 0) {
                return "";
            }
            StringBuilder b = new StringBuilder(str.Length * repeats);
            for (int i = 0; i < repeats; i++) {
                b.Append(str);
            }
            return b.ToString();
        }

    }
}

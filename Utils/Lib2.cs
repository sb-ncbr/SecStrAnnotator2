using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;

namespace SecStrAnnotator2.Utils
{
    public static class Lib2
    {
        public const ConsoleColor ERROR_COLOR = ConsoleColor.Red;
        public const ConsoleColor WARNING_COLOR = ConsoleColor.Red;
        public const ConsoleColor DEBUG_COLOR = ConsoleColor.Blue;
        public static bool DoWriteDebug = true;

        public static void WriteInColor(ConsoleColor color, String s, params object[] args)
        {
            ConsoleColor orig = Console.ForegroundColor;
            Console.ForegroundColor = color;
            Console.Write(String.Format(s, args));
            Console.ForegroundColor = orig;
        }

        public static void WriteDebug(String s, params object[] args)
        {
            if (DoWriteDebug)
            {
                ConsoleColor orig = Console.ForegroundColor;
                Console.ForegroundColor = DEBUG_COLOR;
                Console.Write(String.Format(s, args));
                Console.ForegroundColor = orig;
            }
        }

        public static void WriteLineDebug(String s, params object[] args)
        {
            if (DoWriteDebug)
            {
                ConsoleColor orig = Console.ForegroundColor;
                Console.ForegroundColor = DEBUG_COLOR;
                Console.WriteLine(String.Format(s, args));
                Console.ForegroundColor = orig;
            }
        }

        public static void WriteError(String s, params object[] args)
        {
            ConsoleColor orig = Console.ForegroundColor;
            Console.ForegroundColor = ERROR_COLOR;
            Console.Error.WriteLine("ERROR: " + String.Format(s, args));
            Console.ForegroundColor = orig;
        }

        public static void WriteErrorAndExit(String s, params object[] args)
        {
            WriteError(s, args);
            Environment.Exit(-1);
        }

        public static void WriteWarning(String s, params object[] args)
        {
            ConsoleColor orig = Console.ForegroundColor;
            Console.ForegroundColor = WARNING_COLOR;
            Console.Error.WriteLine("WARNING: " + String.Format(s, args));
            Console.ForegroundColor = orig;
        }

        public static string Enumerate<T>(this IEnumerable<T> list)
        {
            return string.Join(", ", list);
        }
        public static string Enumerate<T>(this IEnumerable<T> list, string separator)
        {
            return string.Join(separator, list);
        }

        public static string EnumerateMultidict<K, V>(this Dictionary<K, List<V>> multidict)
        {
            return multidict.EnumerateMultidict("\n", "\n\t", "\n\t");
        }
        public static string EnumerateMultidict<K, V>(this Dictionary<K, List<V>> multidict, string keySep, string keyValueSep, string valueSep)
        {
            return multidict.Select(kv => kv.Key + keyValueSep + kv.Value.Enumerate(valueSep)).Enumerate(keySep);
        }

        public static IEnumerable<T> WhereHasValue<T>(this IEnumerable<T?> enumer) where T: struct{
            return enumer.Where(x => x.HasValue).Select(x => x.Value);
        }

        public static void Log(object o)
        {
            WriteLineDebug(o.GetType().ToString() + ": " + o.ToString());
        }

        public static void LogList<T>(string name, IEnumerable<T> list)
        {
            WriteLineDebug($"{name} [{list.Count()}]: {list.Enumerate(" ")}");
        }

        /// <summary>
        /// Gets indices on which the value in the array is greater than the previous value (including index 0).
        /// </summary>
        /// <param name="array"> Array of monotonically increasing values.</param>
        public static List<int> RunStartsInOrderedArray(int[] array)
        {
            List<int> runStarts = new List<int>();
            if (array.Length == 0)
            {
                return runStarts;
            }
            else
            {
                runStarts.Add(0);
                int currentRunValue = array[0];
                for (int i = 1; i < array.Length; i++)
                {
                    int value = array[i];
                    // Lib.WriteLineDebug("    " + value + " <--> " + currentRunValue);
                    if (value == currentRunValue)
                    {
                        // do nothing
                    }
                    else if (value > currentRunValue)
                    {
                        // start new run
                        runStarts.Add(i);
                        currentRunValue = value;
                    }
                    else
                    {
                        // array is not ordered (value < currentRunValue)
                        throw new ArgumentException($"Input array was not ordered (array[{i - 1}] = {currentRunValue}, array[{i}] = {value})");
                    }
                }
                return runStarts;
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static T[] AppendAndCopyToArray<T>(this List<T> list, T lastValue)
        {
            int n = list.Count;
            T[] array = new T[n + 1];
            list.CopyTo(array);
            array[n] = lastValue;
            return array;
        }

        public static List<T> Pop<T>(this List<T> queue, int count)
        {
            if (queue.Count < count)
            {
                throw new ArgumentException("Cannot pop " + count.ToString() + " elements from " + queue.Count.ToString());
            }
            List<T> result = queue.GetRange(0, count).ToList();
            queue.RemoveRange(0, count);
            return result;
        }

        public static T Pop<T>(this List<T> queue)
        {
            return queue.Pop(1).First();
        }
    }
}
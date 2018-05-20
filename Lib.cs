using System;
using System.Collections.Generic;
using System.Linq;

namespace SecStrAnnot2
{
    public static class Lib
    {
		public const ConsoleColor ERROR_COLOR = ConsoleColor.Red;
		public const ConsoleColor WARNING_COLOR = ConsoleColor.Red;
		public const ConsoleColor DEBUG_COLOR = ConsoleColor.Blue;
		public static bool DoWriteDebug = true;

		public static void WriteInColor(ConsoleColor color, String s, params object[] args){
			ConsoleColor orig = Console.ForegroundColor;
			Console.ForegroundColor = color;
			Console.Write (String.Format (s,args));
			Console.ForegroundColor = orig;
		}

		public static void WriteDebug(String s, params object[] args){
			if (DoWriteDebug) {
				ConsoleColor orig = Console.ForegroundColor;
				Console.ForegroundColor = DEBUG_COLOR;
				Console.Write (String.Format (s, args));
				Console.ForegroundColor = orig;
			}
		}

		public static void WriteLineDebug(String s, params object[] args){
			if (DoWriteDebug) {
				ConsoleColor orig = Console.ForegroundColor;
				Console.ForegroundColor = DEBUG_COLOR;
				Console.WriteLine (String.Format (s, args));
				Console.ForegroundColor = orig;
			}
		}

		public static void WriteError(String s, params object[] args){
			ConsoleColor orig = Console.ForegroundColor;
			Console.ForegroundColor = ERROR_COLOR;
			Console.Error.WriteLine ("ERROR: " + String.Format (s,args));
			Console.ForegroundColor = orig;
		}

		public static void WriteErrorAndExit(String s, params object[] args){
			WriteError (s, args);
			Environment.Exit (-1);
		}

		public static void WriteWarning(String s, params object[] args){
			ConsoleColor orig = Console.ForegroundColor;
			Console.ForegroundColor = WARNING_COLOR;
			Console.Error.WriteLine ("WARNING: " + String.Format (s,args));
			Console.ForegroundColor = orig;
		}

        public static string Enumerate<T>(this IEnumerable<T> list){
            return string.Join(", ", list);
		}
        public static string Enumerate<T>(this IEnumerable<T> list, string separator){
            return string.Join(separator, list);
		}

        public static string EnumerateMultidict<K,V>(this Dictionary<K,List<V>> multidict){
            return multidict.EnumerateMultidict("\n", "\n\t", "\n\t");
        }
        public static string EnumerateMultidict<K,V>(this Dictionary<K,List<V>> multidict, string keySep, string keyValueSep, string valueSep){
            return multidict.Select(kv => kv.Key + keyValueSep + kv.Value.Enumerate(valueSep)).Enumerate(keySep);
        }

        public static void Log(object o){
            WriteLineDebug(o.GetType().ToString() + ": " + o.ToString());
        }
    }
}
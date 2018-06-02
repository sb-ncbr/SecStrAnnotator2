using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Runtime.Serialization;
using System.Text;
using System.IO;

namespace protein
{
	public enum JsonType{Null,String,Double,Int,Bool,Object,List}
	public class JsonTypeException:Exception{
		public JsonTypeException(String msg) : base(msg) {}
	}
	public class JsonKeyNotFoundException:Exception{
		public String MissingKey { get; private set; }
		public JsonKeyNotFoundException(String msg,String missingKey) : base(msg) {
			MissingKey = missingKey;
		}
	}
	public class JsonIndexException:Exception{
		public int Index { get; private set; }
		public int Count { get; private set; }
		public JsonIndexException(String msg,int index, int count) : base(msg) {
			Index = index;
			Count = count;
		}
	}

	public class JsonValue:IEnumerable<JsonValue>{
		public static String IndentString = "    ";
		public object Content { get; private set; }
		public JsonType Type {
			get {
				if (Content == null)
					return JsonType.Null;
				if (Content is String)
					return JsonType.String;
				if (Content is double)
					return JsonType.Double;
				if (Content is int)
					return JsonType.Int;
				if (Content is Dictionary<string,object>)
					return JsonType.Object;
				if (Content is List<object>)
					return JsonType.List;
				throw new JsonTypeException ("Unknown type: " + Content.GetType ());
			}
		}

		private JsonValue(object content){
			if (content == null || content is String || content is double || content is int || content is bool || content is Dictionary<string,object> || content is List<object>)
				Content = content;
			else
				throw new ArgumentException ("Wrong type of the argument.");
		}

		public JsonValue(String content){
			Content = content;
		}
		public JsonValue(double content){
			if (Double.IsNaN (content) || Double.IsInfinity (content))
				Content = null;
			else
				Content = content;
		}
		public JsonValue(int content){
			Content = content;
		}
		public JsonValue(bool content){
			Content = content;
		}
		public JsonValue(Dictionary<string,object> content){
			Content = content;
		}
		public JsonValue(List<object> content){
			Content = content;
		}
		public JsonValue(){
			Content = null;
		}

		public String String{ 
			get {
				if (Type == JsonType.String)
					return Content as String;
				else
					throw new JsonTypeException ("Not a string.");
			}
		}
		public double Double {
			get {
				if (Type == JsonType.Double)
					return (double)Content;
				if (Type == JsonType.Int)
					return (int)Content;
				else
					throw new JsonTypeException ("Not a double.");
			}
		}
		public int Int {
			get {
				if (Type == JsonType.Int)
					return (int)Content;
				else
					throw new JsonTypeException ("Not an integer.");
			}
		}
		public bool Bool {
			get {
				if (Type == JsonType.Bool)
					return (bool)Content;
				else
					throw new JsonTypeException ("Not a boolean.");
			}
		}
		private Dictionary<string,object> Object {
			get {
				if (Type == JsonType.Object)
					return Content as Dictionary<string,object>;
				else
					throw new JsonTypeException ("Not an object.");
			}
		}
		private List<object> List {
			get {
				if (Type == JsonType.List)
					return Content as List<object>;
				else
					throw new JsonTypeException ("Not a list.");
			}
		}
		public JsonValue this [String key] { 
			get { 
				if (Type == JsonType.Object)
					try {
						return new JsonValue (this.Object [key]);
					} catch (KeyNotFoundException) {
						throw new JsonKeyNotFoundException ("Key \"" + key + "\" not found.", key);
					}
				else
					throw new JsonTypeException ("Not an object.");
			}
			set {
				if (Type == JsonType.Object)
					this.Object [key] = value.Content;
				else
					throw new JsonTypeException ("Not an object.");
			} 
		}
		public void Add(JsonValue value){
			if (Type == JsonType.List)
				(Content as List<object>).Add (value.Content);
			else
				throw new JsonTypeException ("Not a list.");
		}
		public JsonValue this [int i] { 
			get { 
				if (Type == JsonType.List)
					try {
						return new JsonValue (this.List [i]);
					} catch (ArgumentOutOfRangeException) {
						throw new JsonIndexException ("Index out of range (" + i + " of " + this.Count + ").", i, this.Count);
					}
				else if (Type == JsonType.Object)
					try {
						return new JsonValue (this.Object.ElementAt (i).Value);
					} catch (ArgumentOutOfRangeException) {
						throw new JsonIndexException ("Index out of range (" + i + " of " + this.Count + ").", i, this.Count);
					}
				else
					throw new JsonTypeException ("Not a list.");
			}
			set {
				if (Type == JsonType.List)
					this.List [i] = value.Content;
				else
					throw new JsonTypeException ("Not a list.");
			} 
		}
		public bool Contains (string key){
			if (Type == JsonType.Object)
				return this.Object.ContainsKey (key);
			else
				throw new JsonTypeException ("Not an object.");
		}
		public int Count {
			get {
				if (Type == JsonType.List)
					return this.List.Count;
				else if (Type == JsonType.Object)
					return this.Object.Count;
				else
					throw new JsonTypeException ("Not an object.");
			}
		}
		public string Key(int i){
			if (Type == JsonType.Object)
				return this.Object.ElementAt (i).Key;
			else
				throw new JsonTypeException ("Not an object.");
		}

		public IEnumerator<JsonValue> GetEnumerator(){
			if (Type == JsonType.List)
				return new JsonEnumerator (this.List.GetEnumerator ());
			else
				throw new JsonTypeException ("Not a list.");
		}
		IEnumerator System.Collections.IEnumerable.GetEnumerator(){
			return GetEnumerator ();
		}
		public class JsonEnumerator : IEnumerator<JsonValue>{
			private IEnumerator<object> inner;
			public JsonEnumerator(IEnumerator<object> innerEnumerator){
				inner=innerEnumerator;
			}

			public JsonValue Current{ get { return new JsonValue (inner.Current); } }
			object System.Collections.IEnumerator.Current{ get { return Current; } }
			public void Dispose(){ inner.Dispose (); }
			public bool MoveNext(){ return inner.MoveNext (); }
			public void Reset (){ inner.Reset (); }
		}

		public static JsonValue MakeList ()
		{
			return new JsonValue (new List<object>());
		}
		public static JsonValue MakeList (IEnumerable<JsonValue> list)
		{
			List<object> cont = list.Select (j => j.Content).ToList ();
			return new JsonValue (cont);
		}
		public static JsonValue MakeObject (params object[] fields)
		{
			if (fields.Length%2!=0)
				throw new ArgumentException ("Odd number of arguments.");
			Dictionary<string,object> cont = new Dictionary<string,object> ();
			for (int i = 0; i < fields.Length/2; i++) {
				object key = fields [2 * i];
				object value = fields [2 * i + 1];
				if (key is string || value is JsonValue)
					cont [key as string] = (value as JsonValue).Content;
				else
					throw new ArgumentException ("Odd arguments must be strings, even arguments JsonValues.");
			}
			return new JsonValue (cont);
		}

		public static JsonValue FromString(String str){
			object obj = JSONParser.FromJson<object> (str);
			return new JsonValue (obj);
		}
		public static JsonValue FromFile(String filename){
			String content = "{}";
			using (StreamReader r = new StreamReader (filename)) {
				content = r.ReadToEnd ();
			}
			JsonValue value = FromString (content);
			return value;
		}
		public override String ToString(){
			// Use this Firefox add-on to show Json in a pretty way: https://addons.mozilla.org/en-US/firefox/addon/jsonview/
			return JSONPrettyWriter.ToJson (this.Content, 0);
		}
		public String ToString(int maxIndentLevel){
			// Use this Firefox add-on to show Json in a pretty way: https://addons.mozilla.org/en-US/firefox/addon/jsonview/
			return JSONPrettyWriter.ToJson (this.Content, maxIndentLevel);
		}

		// The following 2 nested classes were taken from https://github.com/zanders3/json

		// Really simple JSON parser in ~300 lines
		// - Attempts to parse JSON files with minimal GC allocation
		// - Nice and simple "[1,2,3]".FromJson<List<int>>() API
		// - Classes and structs can be parsed too!
		//      class Foo { public int Value; }
		//      "{\"Value\":10}".FromJson<Foo>()
		// - Can parse JSON without type information into Dictionary<string,object> and List<object> e.g.
		//      "[1,2,3]".FromJson<object>().GetType() == typeof(List<object>)
		//      "{\"Value\":10}".FromJson<object>().GetType() == typeof(Dictionary<string,object>)
		// - No JIT Emit support to support AOT compilation on iOS
		// - Attempts are made to NOT throw an exception if the JSON is corrupted or invalid: returns null instead.
		// - Only public fields and property setters on classes/structs will be written to
		//
		// Limitations:
		// - No JIT Emit support to parse structures quickly
		// - Limited to parsing <2GB JSON files (due to int.MaxValue)
		// - Parsing of abstract classes or interfaces is NOT supported and will throw an exception.
		private static class JSONParser
		{
			static Stack<List<string>> splitArrayPool = new Stack<List<string>>();
			static StringBuilder stringBuilder = new StringBuilder();
			static readonly Dictionary<Type, Dictionary<string, FieldInfo>> fieldInfoCache = new Dictionary<Type, Dictionary<string, FieldInfo>>();
			static readonly Dictionary<Type, Dictionary<string, PropertyInfo>> propertyInfoCache = new Dictionary<Type, Dictionary<string, PropertyInfo>>();

			public static T FromJson<T>(string json)
			{
				//Remove all whitespace not within strings to make parsing simpler
				stringBuilder.Length = 0;
				for (int i = 0; i < json.Length; i++)
				{
					char c = json[i];
					if (c == '\"')
					{
						i = AppendUntilStringEnd(true, i, json);
						continue;
					}
					if (char.IsWhiteSpace(c))
						continue;

					stringBuilder.Append(c);
				}

				//Parse the thing!
				return (T)ParseValue(typeof(T), stringBuilder.ToString());
			}

			static int AppendUntilStringEnd(bool appendEscapeCharacter, int startIdx, string json)
			{
				stringBuilder.Append(json[startIdx]);
				for (int i = startIdx+1; i<json.Length; i++)
				{
					if (json[i] == '\\')
					{
						if (appendEscapeCharacter)
							stringBuilder.Append(json[i]);
						stringBuilder.Append(json[i + 1]);
						i++;//Skip next character as it is escaped
					}
					else if (json[i] == '\"')
					{
						stringBuilder.Append(json[i]);
						return i;
					}
					else
						stringBuilder.Append(json[i]);
				}
				return json.Length - 1;
			}

			//Splits { <value>:<value>, <value>:<value> } and [ <value>, <value> ] into a list of <value> strings
			static List<string> Split(string json)
			{
				List<string> splitArray = splitArrayPool.Count > 0 ? splitArrayPool.Pop() : new List<string>();
				splitArray.Clear();
				int parseDepth = 0;
				stringBuilder.Length = 0;
				for (int i = 1; i<json.Length-1; i++)
				{
					switch (json[i])
					{
					case '[':
					case '{':
						parseDepth++;
						break;
					case ']':
					case '}':
						parseDepth--;
						break;
					case '\"':
						i = AppendUntilStringEnd(true, i, json);
						continue;
					case ',':
					case ':':
						if (parseDepth == 0)
						{
							splitArray.Add(stringBuilder.ToString());
							stringBuilder.Length = 0;
							continue;
						}
						break;
					}

					stringBuilder.Append(json[i]);
				}

				splitArray.Add(stringBuilder.ToString());

				return splitArray;
			}

			internal static object ParseValue(Type type, string json)
			{
				if (type == typeof(string))
				{
					if (json.Length <= 2)
						return string.Empty;
					string str = json.Substring(1, json.Length - 2);
					return str.Replace("\\\\", "\"\"").Replace("\\", string.Empty).Replace("\"\"", "\\");
				}
				if (type == typeof(int))
				{
					int result;
					int.TryParse(json, out result);
					return result;
				}
				if (type == typeof(float))
				{
					float result;
					float.TryParse(json, out result);
					return result;
				}
				if (type == typeof(double))
				{
					double result;
					double.TryParse(json, out result);
					return result;
				}
				if (type == typeof(bool))
				{
					return json.ToLower() == "true";
				}
				if (json == "null")
				{
					return null;
				}
				if (type.IsArray)
				{
					Type arrayType = type.GetElementType();
					if (json[0] != '[' || json[json.Length - 1] != ']')
						return null;

					List<string> elems = Split(json);
					Array newArray = Array.CreateInstance(arrayType, elems.Count);
					for (int i = 0; i < elems.Count; i++)
						newArray.SetValue(ParseValue(arrayType, elems[i]), i);
					splitArrayPool.Push(elems);
					return newArray;
				}
				if (type.IsGenericType && type.GetGenericTypeDefinition() == typeof(List<>))
				{
					Type listType = type.GetGenericArguments()[0];
					if (json[0] != '[' || json[json.Length - 1] != ']')
						return null;

					List<string> elems = Split(json);
					var list = (IList)type.GetConstructor(new Type[] { typeof(int) }).Invoke(new object[] { elems.Count });
					for (int i = 0; i < elems.Count; i++)
						list.Add(ParseValue(listType, elems[i]));
					splitArrayPool.Push(elems);
					return list;
				}
				if (type.IsGenericType && type.GetGenericTypeDefinition() == typeof(Dictionary<,>))
				{
					Type keyType, valueType;
					{
						Type[] args = type.GetGenericArguments();
						keyType = args[0];
						valueType = args[1];
					}

					//Refuse to parse dictionary keys that aren't of type string
					if (keyType != typeof(string))
						return null;
					//Must be a valid dictionary element
					if (json[0] != '{' || json[json.Length - 1] != '}')
						return null;
					//The list is split into key/value pairs only, this means the split must be divisible by 2 to be valid JSON
					List<string> elems = Split(json);
					if (elems.Count % 2 != 0)
						return null;

					var dictionary = (IDictionary)type.GetConstructor(new Type[] { typeof(int) }).Invoke(new object[] { elems.Count / 2 });
					for (int i = 0; i < elems.Count; i += 2)
					{
						if (elems[i].Length <= 2)
							continue;
						string keyValue = elems[i].Substring(1, elems[i].Length - 2);
						object val = ParseValue(valueType, elems[i + 1]);
						dictionary.Add(keyValue, val);
					}
					return dictionary;
				}
				if (type == typeof(object))
				{
					return ParseAnonymousValue(json);
				}    
				if (json[0] == '{' && json[json.Length - 1] == '}')
				{
					return ParseObject(type, json);
				}

				return null;
			}

			static object ParseAnonymousValue(string json)
			{
				if (json.Length == 0)
					return null;
				if (json[0] == '{' && json[json.Length - 1] == '}')
				{
					List<string> elems = Split(json);
					if (elems.Count % 2 != 0)
						return null;
					var dict = new Dictionary<string, object>(elems.Count / 2);
					for (int i = 0; i < elems.Count; i += 2)
						dict.Add(elems[i].Substring(1, elems[i].Length - 2), ParseAnonymousValue(elems[i + 1]));
					return dict;
				}
				if (json[0] == '[' && json[json.Length - 1] == ']')
				{
					List<string> items = Split(json);
					var finalList = new List<object>(items.Count);
					//Console.WriteLine ("JSON parsing: " + items.Count.ToString () + ": '" + items.EnumerateWithCommas () + "'");
					if (items.Count == 1 && items[0]=="") { //@Adam: to avoid empty array [] being read as [null]
						return finalList;
					}
					for (int i = 0; i < items.Count; i++)
						finalList.Add(ParseAnonymousValue(items[i]));
					return finalList;
				}
				if (json[0] == '\"' && json[json.Length - 1] == '\"')
				{
					string str = json.Substring(1, json.Length - 2);
					return str.Replace("\\\\", "\"\"").Replace("\\", string.Empty).Replace("\"\"", "\\");
				}
				if (char.IsDigit(json[0]) || json[0] == '-')
				{
					if (json.Contains("."))
					{
						double result;
						double.TryParse(json, out result);
						return result;
					}
					else
					{
						int result;
						int.TryParse(json, out result);
						return result;
					}
				}
				if (json == "true")
					return true;
				if (json == "false")
					return false;
				// handles json == "null" as well as invalid JSON
				return null;
			}

			static object ParseObject(Type type, string json)
			{
				object instance = FormatterServices.GetUninitializedObject(type);

				//The list is split into key/value pairs only, this means the split must be divisible by 2 to be valid JSON
				List<string> elems = Split(json);
				if (elems.Count % 2 != 0)
					return instance;

				Dictionary<string, FieldInfo> nameToField;
				Dictionary<string, PropertyInfo> nameToProperty;
				if (!fieldInfoCache.TryGetValue(type, out nameToField))
				{
					nameToField = type.GetFields().Where(field => field.IsPublic).ToDictionary(field => field.Name);
					fieldInfoCache.Add(type, nameToField);
				}
				if (!propertyInfoCache.TryGetValue(type, out nameToProperty))
				{
					nameToProperty = type.GetProperties().ToDictionary(p => p.Name);
					propertyInfoCache.Add(type, nameToProperty);
				}

				for (int i = 0; i < elems.Count; i += 2)
				{
					if (elems[i].Length <= 2)
						continue;
					string key = elems[i].Substring(1, elems[i].Length - 2);
					string value = elems[i + 1];

					FieldInfo fieldInfo;
					PropertyInfo propertyInfo;
					if (nameToField.TryGetValue(key, out fieldInfo))
						fieldInfo.SetValue(instance, ParseValue(fieldInfo.FieldType, value));
					else if (nameToProperty.TryGetValue(key, out propertyInfo))
						propertyInfo.SetValue(instance, ParseValue(propertyInfo.PropertyType, value), null);
				}

				return instance;
			}
		}

		//Really simple JSON writer
		//- Outputs JSON structures from an object
		//- Really simple API (new List<int> { 1, 2, 3 }).ToJson() == "[1,2,3]"
		//- Will only output public fields and property getters on objects
		private static class JSONWriter
		{
			public static string ToJson(object item)
			{
				StringBuilder stringBuilder = new StringBuilder();
				AppendValue(stringBuilder, item);
				return stringBuilder.ToString();
			}

			static void AppendValue(StringBuilder stringBuilder, object item)
			{
				if (item == null)
				{
					stringBuilder.Append("null");
					return;
				}

				Type type = item.GetType();
				if (type == typeof(string))
				{
					stringBuilder.Append('\"');
					stringBuilder.Append(((string)item).Replace("\\", "\\\\"));
					stringBuilder.Append('\"');
				}
				else if (type == typeof(int) || type == typeof(float) || type == typeof(double))
				{
					stringBuilder.Append(item.ToString());
				}
				else if (type == typeof(bool))
				{
					stringBuilder.Append(((bool)item) ? "true" : "false");
				}
				else if (item is IList)
				{
					stringBuilder.Append('[');
					bool isFirst = true;
					IList list = item as IList;
					for (int i = 0; i < list.Count; i++)
					{
						if (isFirst)
							isFirst = false;
						else
							stringBuilder.Append(',');
						AppendValue(stringBuilder, list[i]);
					}
					stringBuilder.Append(']');
				}
				else if (type.IsGenericType && type.GetGenericTypeDefinition() == typeof(Dictionary<,>))
				{
					Type keyType = type.GetGenericArguments()[0];

					//Refuse to output dictionary keys that aren't of type string
					if (keyType != typeof(string))
					{
						stringBuilder.Append("{}");
						return;
					}

					stringBuilder.Append('{');
					IDictionary dict = item as IDictionary;
					bool isFirst = true;
					foreach (object key in dict.Keys)
					{
						if (isFirst)
							isFirst = false;
						else
							stringBuilder.Append(',');
						stringBuilder.Append('\"');
						stringBuilder.Append((string)key);
						stringBuilder.Append("\":");
						AppendValue(stringBuilder, dict[key]);
					}
					stringBuilder.Append('}');
				}
				else
				{
					stringBuilder.Append('{');

					bool isFirst = true;
					FieldInfo[] fieldInfos = type.GetFields();
					for (int i = 0; i < fieldInfos.Length; i++)
					{
						if (fieldInfos[i].IsPublic)
						{
							object value = fieldInfos[i].GetValue(item);
							if (value != null)
							{
								if (isFirst)
									isFirst = false;
								else
									stringBuilder.Append(',');
								stringBuilder.Append('\"');
								stringBuilder.Append(fieldInfos[i].Name);
								stringBuilder.Append("\":");
								AppendValue(stringBuilder, value);
							}
						}
					}
					PropertyInfo[] propertyInfo = type.GetProperties();
					for (int i = 0; i<propertyInfo.Length; i++)
					{
						if (propertyInfo[i].CanRead)
						{
							object value = propertyInfo[i].GetValue(item, null);
							if (value != null)
							{
								if (isFirst)
									isFirst = false;
								else
									stringBuilder.Append(',');
								stringBuilder.Append('\"');
								stringBuilder.Append(propertyInfo[i].Name);
								stringBuilder.Append("\":");
								AppendValue(stringBuilder, value);
							}
						}
					}

					stringBuilder.Append('}');
				}
			}
		}


		//Really simple JSON writer
		//- Outputs JSON structures from an object
		//- Really simple API (new List<int> { 1, 2, 3 }).ToJson() == "[1,2,3]"
		//- Will only output public fields and property getters on objects
		private static class JSONPrettyWriter
		{
			
			public static string ToJson(object item, int maxIndentLevel)
			{
				StringBuilder stringBuilder = new StringBuilder();
				AppendValue(stringBuilder, item, 0, maxIndentLevel);
				return stringBuilder.ToString();
			}

			private static void AppendIndent(StringBuilder stringBuilder, int indent, int maxIndent){
				// Print indent
				if (indent <= maxIndent) {
					stringBuilder.Append ('\n');
					for (int i = 0; i < indent; i++) {
						stringBuilder.Append (JsonValue.IndentString);
					}
				} else {
					stringBuilder.Append (' ');
				}
			}

			static void AppendValue(StringBuilder stringBuilder, object item, int currentIndent, int maxIndent)
			{
				// Print value
				if (item == null)
				{
					// Null
					stringBuilder.Append("null");
					return;
				}

				Type type = item.GetType();
				if (type == typeof(string))
				{
					// String
					stringBuilder.Append('\"');
					stringBuilder.Append(((string)item).Replace("\\", "\\\\"));
					stringBuilder.Append('\"');
				}
				else if (type == typeof(int) || type == typeof(float) || type == typeof(double))
				{
					// Number
					stringBuilder.Append(item.ToString());
				}
				else if (type == typeof(bool))
				{
					// True, False
					stringBuilder.Append(((bool)item) ? "true" : "false");
				}
				else if (item is IList)
				{
					// Array
					IList list = item as IList;
					if (list.Count == 0) {
						// Empty array
						stringBuilder.Append ("[ ]");
					} else {
						// Non-empty array
						stringBuilder.Append ('[');
						bool isFirst = true;
						for (int i = 0; i < list.Count; i++) {
							if (i > 0) {
								stringBuilder.Append (',');
							}
							AppendIndent(stringBuilder, currentIndent + 1, maxIndent);
							AppendValue (stringBuilder, list [i], currentIndent + 1, maxIndent);
						}
						if (currentIndent < maxIndent) {
							AppendIndent (stringBuilder, currentIndent, maxIndent);
						} else {
							stringBuilder.Append (' ');
						}
						stringBuilder.Append (']');
					}
				}
				else if (type.IsGenericType && type.GetGenericTypeDefinition() == typeof(Dictionary<,>))
				{
					// Object (dictionary)
					Type keyType = type.GetGenericArguments()[0];

					//Refuse to output dictionary keys that aren't of type string
					if (keyType != typeof(string))
					{
						stringBuilder.Append("{ }");
						return;
					}

					IDictionary dict = item as IDictionary;
					if (dict.Count == 0) {
						// Empty dictionary
						stringBuilder.Append ("{ }");
					} else {
						// Non-empty dictionary
						stringBuilder.Append ('{');
						bool isFirst = true;
						foreach (object key in dict.Keys) {
							if (isFirst)
								isFirst = false;
							else
								stringBuilder.Append (',');
							AppendIndent(stringBuilder, currentIndent + 1, maxIndent);
							stringBuilder.Append ('\"');
							stringBuilder.Append ((string)key);
							stringBuilder.Append ("\": ");
							AppendValue (stringBuilder, dict [key], currentIndent + 1, maxIndent);
						}
						if (currentIndent < maxIndent) {
							AppendIndent (stringBuilder, currentIndent, maxIndent);
						} else {
							stringBuilder.Append (' ');
						}
						stringBuilder.Append ('}');
					}
				}
				else
				{
					// Object (object)
					stringBuilder.Append('{');

					bool isFirst = true;
					FieldInfo[] fieldInfos = type.GetFields();
					for (int i = 0; i < fieldInfos.Length; i++)
					{
						if (fieldInfos[i].IsPublic)
						{
							object value = fieldInfos[i].GetValue(item);
							if (value != null)
							{
								if (isFirst)
									isFirst = false;
								else
									stringBuilder.Append(',');
								AppendValue (stringBuilder, value, currentIndent + 1, maxIndent);
								stringBuilder.Append('\"');
								stringBuilder.Append(fieldInfos[i].Name);
								stringBuilder.Append("\":");
								AppendValue(stringBuilder, value, currentIndent + 1, maxIndent);
							}
						}
					}
					PropertyInfo[] propertyInfo = type.GetProperties();
					for (int i = 0; i<propertyInfo.Length; i++)
					{
						if (propertyInfo[i].CanRead)
						{
							object value = propertyInfo[i].GetValue(item, null);
							if (value != null)
							{
								if (isFirst)
									isFirst = false;
								else
									stringBuilder.Append(',');
								AppendValue (stringBuilder, value, currentIndent + 1, maxIndent);
								stringBuilder.Append('\"');
								stringBuilder.Append(propertyInfo[i].Name);
								stringBuilder.Append("\":");
								AppendValue(stringBuilder, value, currentIndent + 1, maxIndent);
							}
						}
					}
					if (currentIndent < maxIndent) {
						AppendIndent (stringBuilder, currentIndent, maxIndent);
					} else {
						stringBuilder.Append (' ');
					}
					stringBuilder.Append('}');
				}
				if (currentIndent == 0 && maxIndent > 0) {
					stringBuilder.Append ('\n');
				}
			}
		}

	}


	}
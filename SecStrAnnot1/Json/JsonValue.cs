using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Runtime.Serialization;
using System.Text;
using System.IO;

using protein.Json.Helpers;

namespace protein.Json
{
    public class JsonValue : IEnumerable<JsonValue>
    {
        public static String IndentString = "    ";

        public object Content { get; private set; }
        
        public JsonType Type
        {
            get
            {
                if (Content == null)
                    return JsonType.Null;
                if (Content is String)
                    return JsonType.String;
                if (Content is double)
                    return JsonType.Double;
                if (Content is int)
                    return JsonType.Int;
                if (Content is Dictionary<string, object>)
                    return JsonType.Object;
                if (Content is List<object>)
                    return JsonType.List;
                throw new JsonTypeException("Unknown type: " + Content.GetType());
            }
        }

        internal JsonValue(object content)
        {
            if (content == null || content is String || content is double || content is int || content is bool || content is Dictionary<string, object> || content is List<object>)
                Content = content;
            else
                throw new ArgumentException("Wrong type of the argument.");
        }

        public JsonValue(String content)
        {
            Content = content;
        }

        public JsonValue(double content)
        {
            if (Double.IsNaN(content) || Double.IsInfinity(content))
                Content = null;
            else
                Content = content;
        }
        
        public JsonValue(int content)
        {
            Content = content;
        }
        
        public JsonValue(bool content)
        {
            Content = content;
        }
        
        public JsonValue(Dictionary<string, object> content)
        {
            Content = content;
        }
        
        public JsonValue(List<object> content)
        {
            Content = content;
        }
        
        public JsonValue()
        {
            Content = null;
        }

        public String String
        {
            get
            {
                if (Type == JsonType.String)
                    return Content as String;
                else
                    throw new JsonTypeException("Not a string.");
            }
        }

        public double Double
        {
            get
            {
                if (Type == JsonType.Double)
                    return (double)Content;
                if (Type == JsonType.Int)
                    return (int)Content;
                else
                    throw new JsonTypeException("Not a double.");
            }
        }

        public int Int
        {
            get
            {
                if (Type == JsonType.Int)
                    return (int)Content;
                else
                    throw new JsonTypeException("Not an integer.");
            }
        }

        public bool Bool
        {
            get
            {
                if (Type == JsonType.Bool)
                    return (bool)Content;
                else
                    throw new JsonTypeException("Not a boolean.");
            }
        }

        private Dictionary<string, object> Object
        {
            get
            {
                if (Type == JsonType.Object)
                    return Content as Dictionary<string, object>;
                else
                    throw new JsonTypeException("Not an object.");
            }
        }

        private List<object> List
        {
            get
            {
                if (Type == JsonType.List)
                    return Content as List<object>;
                else
                    throw new JsonTypeException("Not a list.");
            }
        }

        public JsonValue this[String key]
        {
            get
            {
                if (Type == JsonType.Object)
                    try
                    {
                        return new JsonValue(this.Object[key]);
                    }
                    catch (KeyNotFoundException)
                    {
                        throw new JsonKeyNotFoundException("Key \"" + key + "\" not found.", key);
                    }
                else
                    throw new JsonTypeException("Not an object.");
            }
            set
            {
                if (Type == JsonType.Object)
                    this.Object[key] = value.Content;
                else
                    throw new JsonTypeException("Not an object.");
            }
        }

        /**Append a value at the end of the JSON list.*/
        public void Add(JsonValue value)
        {
            if (Type == JsonType.List)
                (Content as List<object>).Add(value.Content);
            else
                throw new JsonTypeException("Not a list.");
        }

        public JsonValue this[int i]
        {
            get
            {
                if (Type == JsonType.List)
                    try
                    {
                        return new JsonValue(this.List[i]);
                    }
                    catch (ArgumentOutOfRangeException)
                    {
                        throw new JsonIndexException("Index out of range (" + i + " of " + this.Count + ").", i, this.Count);
                    }
                else if (Type == JsonType.Object)
                    try
                    {
                        return new JsonValue(this.Object.ElementAt(i).Value);
                    }
                    catch (ArgumentOutOfRangeException)
                    {
                        throw new JsonIndexException("Index out of range (" + i + " of " + this.Count + ").", i, this.Count);
                    }
                else
                    throw new JsonTypeException("Not a list.");
            }
            set
            {
                if (Type == JsonType.List)
                    this.List[i] = value.Content;
                else
                    throw new JsonTypeException("Not a list.");
            }
        }

        public bool Contains(string key)
        {
            if (Type == JsonType.Object)
                return this.Object.ContainsKey(key);
            else
                throw new JsonTypeException("Not an object.");
        }

        public int Count
        {
            get
            {
                if (Type == JsonType.List)
                    return this.List.Count;
                else if (Type == JsonType.Object)
                    return this.Object.Count;
                else
                    throw new JsonTypeException("Not an object.");
            }
        }

        public string Key(int i)
        {
            if (Type == JsonType.Object)
                return this.Object.ElementAt(i).Key;
            else
                throw new JsonTypeException("Not an object.");
        }

        public IEnumerator<JsonValue> GetEnumerator()
        {
            if (Type == JsonType.List)
                return new JsonEnumerator(this.List.GetEnumerator());
            else
                throw new JsonTypeException("Not a list.");
        }

        IEnumerator System.Collections.IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }


        public static JsonValue MakeList()
        {
            return new JsonValue(new List<object>());
        }

        public static JsonValue MakeList(IEnumerable<JsonValue> list)
        {
            List<object> cont = list.Select(j => j.Content).ToList();
            return new JsonValue(cont);
        }

        public static JsonValue MakeObject(params object[] fields)
        {
            if (fields.Length % 2 != 0)
                throw new ArgumentException("Odd number of arguments.");
            Dictionary<string, object> cont = new Dictionary<string, object>();
            for (int i = 0; i < fields.Length / 2; i++)
            {
                object key = fields[2 * i];
                object value = fields[2 * i + 1];
                if (key is string || value is JsonValue)
                    cont[key as string] = (value as JsonValue).Content;
                else
                    throw new ArgumentException("Odd arguments must be strings, even arguments JsonValues.");
            }
            return new JsonValue(cont);
        }

        public static JsonValue FromString(String str)
        {
            object obj = JSONParser.FromJson<object>(str);
            return new JsonValue(obj);
        }

        public static JsonValue FromFile(String filename)
        {
            String content = "{}";
            using (StreamReader r = new StreamReader(filename))
            {
                content = r.ReadToEnd();
            }
            JsonValue value = FromString(content);
            return value;
        }

        public override String ToString()
        {
            // Use this Firefox add-on to show Json in a pretty way: https://addons.mozilla.org/en-US/firefox/addon/jsonview/
            return JsonPrettyWriter.ToJson(this.Content, 0);
        }

        public String ToString(int maxIndentLevel)
        {
            // Use this Firefox add-on to show Json in a pretty way: https://addons.mozilla.org/en-US/firefox/addon/jsonview/
            return JsonPrettyWriter.ToJson(this.Content, maxIndentLevel);
        }

    }
}
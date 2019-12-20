using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Runtime.Serialization;
using System.Text;
using System.IO;

namespace protein.Json.Helpers
{
    // Taken from https://github.com/zanders3/json

    //Really simple JSON writer
    //- Outputs JSON structures from an object
    //- Really simple API (new List<int> { 1, 2, 3 }).ToJson() == "[1,2,3]"
    //- Will only output public fields and property getters on objects

    static class JsonPrettyWriter
    {

        public static string ToJson(object item, int maxIndentLevel)
        {
            StringBuilder stringBuilder = new StringBuilder();
            AppendValue(stringBuilder, item, 0, maxIndentLevel);
            return stringBuilder.ToString();
        }

        private static void AppendIndent(StringBuilder stringBuilder, int indent, int maxIndent)
        {
            // Print indent
            if (indent <= maxIndent)
            {
                stringBuilder.Append('\n');
                for (int i = 0; i < indent; i++)
                {
                    stringBuilder.Append(JsonValue.IndentString);
                }
            }
            else
            {
                stringBuilder.Append(' ');
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
                if (list.Count == 0)
                {
                    // Empty array
                    stringBuilder.Append("[ ]");
                }
                else
                {
                    // Non-empty array
                    stringBuilder.Append('[');
                    // bool isFirst = true;
                    for (int i = 0; i < list.Count; i++)
                    {
                        if (i > 0)
                        {
                            stringBuilder.Append(',');
                        }
                        AppendIndent(stringBuilder, currentIndent + 1, maxIndent);
                        AppendValue(stringBuilder, list[i], currentIndent + 1, maxIndent);
                    }
                    if (currentIndent < maxIndent)
                    {
                        AppendIndent(stringBuilder, currentIndent, maxIndent);
                    }
                    else
                    {
                        stringBuilder.Append(' ');
                    }
                    stringBuilder.Append(']');
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
                if (dict.Count == 0)
                {
                    // Empty dictionary
                    stringBuilder.Append("{ }");
                }
                else
                {
                    // Non-empty dictionary
                    stringBuilder.Append('{');
                    bool isFirst = true;
                    foreach (object key in dict.Keys)
                    {
                        if (isFirst)
                            isFirst = false;
                        else
                            stringBuilder.Append(',');
                        AppendIndent(stringBuilder, currentIndent + 1, maxIndent);
                        stringBuilder.Append('\"');
                        stringBuilder.Append((string)key);
                        stringBuilder.Append("\": ");
                        AppendValue(stringBuilder, dict[key], currentIndent + 1, maxIndent);
                    }
                    if (currentIndent < maxIndent)
                    {
                        AppendIndent(stringBuilder, currentIndent, maxIndent);
                    }
                    else
                    {
                        stringBuilder.Append(' ');
                    }
                    stringBuilder.Append('}');
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
                            AppendValue(stringBuilder, value, currentIndent + 1, maxIndent);
                            stringBuilder.Append('\"');
                            stringBuilder.Append(fieldInfos[i].Name);
                            stringBuilder.Append("\":");
                            AppendValue(stringBuilder, value, currentIndent + 1, maxIndent);
                        }
                    }
                }
                PropertyInfo[] propertyInfo = type.GetProperties();
                for (int i = 0; i < propertyInfo.Length; i++)
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
                            AppendValue(stringBuilder, value, currentIndent + 1, maxIndent);
                            stringBuilder.Append('\"');
                            stringBuilder.Append(propertyInfo[i].Name);
                            stringBuilder.Append("\":");
                            AppendValue(stringBuilder, value, currentIndent + 1, maxIndent);
                        }
                    }
                }
                if (currentIndent < maxIndent)
                {
                    AppendIndent(stringBuilder, currentIndent, maxIndent);
                }
                else
                {
                    stringBuilder.Append(' ');
                }
                stringBuilder.Append('}');
            }
            if (currentIndent == 0 && maxIndent > 0)
            {
                stringBuilder.Append('\n');
            }
        }

    }
}
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace logger
{
    class Logger
    {
        TextWriter writer;
        List<String> blocks;
        public Level LogLevel { get; set; }
        public bool TotalIgnore { get; set; }

        public Logger(TextWriter writer, bool totalIgnore)
        {
            if (writer == null) throw new ArgumentNullException();
            this.writer = writer;
            this.blocks = new List<String>();
            this.LogLevel = Level.INFO;
            this.TotalIgnore = totalIgnore;
            if (!TotalIgnore)
            {
                writer.WriteLine("<?xml version=\"1.0\" encoding=\"utf-8\"?>");
                OpenBlock("log");
            }
        }

        public Logger(TextWriter writer) : this(writer, false) { }

        public Logger(TextWriter writer, Level level)
            : this(writer)
        {
            this.LogLevel = level;
        }

        public void Close()
        {
            if (blocks.Count != 1) throw new InvalidOperationException("Not all blocks are closed.");
            CloseBlock("log");
            writer.Close();
            writer = null;
        }

        public void CloseFatal()
        {
            writer.Close();
            writer = null;
        }

        public void Debug(String msg)
        {
            if (LogLevel >= Level.DEBUG && !TotalIgnore)
            {
                for (int i = 0; i < blocks.Count; i++)
                {
                    writer.Write("\t");
                }
                writer.WriteLine("<debug>" + msg + "</debug>");
            }
        }

        public void Info(String msg)
        {
            if (this.LogLevel >= Level.INFO && !TotalIgnore)
            {
                for (int i = 0; i < blocks.Count; i++)
                {
                    writer.Write("\t");
                }
                writer.WriteLine("<info>" + msg + "</info>");
            }
        }
        
        public void Error(String msg)
        {
            if (this.LogLevel >= Level.ERROR && !TotalIgnore)
            {
                for (int i = 0; i < blocks.Count; i++)
                {
                    writer.Write("\t");
                }
                writer.WriteLine("<error>" + msg + "</error>");
            }
        }

        public void Fatal(String msg)
        {
            if (LogLevel >= Level.FATAL && !TotalIgnore)
            {
                writer.WriteLine("<fatal>" + msg + "</fatal>");
            }
        }

        public void OpenBlockD(String blockName)
        {
            if (this.LogLevel >= Level.DEBUG && !TotalIgnore)
            {
                for (int i = 0; i < blocks.Count; i++)
                {
                    writer.Write("\t");
                }
                writer.WriteLine("<blockD name=\"" + blockName + "\">");
                blocks.Add(blockName);
            }
        }

        public void CloseBlockD(String blockName)
        {
            if (this.LogLevel >= Level.DEBUG && !TotalIgnore)
            {
                if (blockName.Equals(blocks[blocks.Count - 1]))
                {
                    blocks.RemoveAt(blocks.Count - 1);
                    for (int i = 0; i < blocks.Count; i++)
                    {
                        writer.Write("\t");
                    }
                    writer.WriteLine("</blockD>");
                }
                else
                {
                    throw new InvalidOperationException("Opened block \"" + blocks[blocks.Count - 1] + "\"" + "closing block \"" + blockName + "\"");
                }
            }
        }
       
        public void OpenBlock(String blockName)
        {
            if (this.LogLevel >= Level.INFO && !TotalIgnore)
            {
                for (int i = 0; i < blocks.Count; i++)
                {
                    writer.Write("\t");
                }
                writer.WriteLine("<block name=\"" + blockName + "\">");
                blocks.Add(blockName);
            }
        }

        public void CloseBlock(String blockName)
        {
            if (this.LogLevel >= Level.INFO && !TotalIgnore)
            {
                if (blockName.Equals(blocks[blocks.Count - 1]))
                {
                    blocks.RemoveAt(blocks.Count - 1);
                    for (int i = 0; i < blocks.Count; i++)
                    {
                        writer.Write("\t");
                    }
                    writer.WriteLine("</block>");
                }
                else
                {
                    throw new InvalidOperationException("Opened block \"" + blocks[blocks.Count - 1] + "\"" + "closing block \"" + blockName + "\"");
                }
            }
        }

        public enum Level { FATAL=1, ERROR, INFO, DEBUG}
    }
}

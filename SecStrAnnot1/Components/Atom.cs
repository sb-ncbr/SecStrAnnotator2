using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Globalization;

using protein.Geometry;

namespace protein.Components
{
    public class Atom : IComparable<Atom>
    {
        //http://www.wwpdb.org/documentation/format33/sect9.html#ATOM

        public int Serial { get; private set; }

        private String name;
        public String Name
        {
            get
            {
                return name;
            }
            private set
            {
                name = value;
                // if (value.Length == 4)
                // {
                //     name = value;
                // }
                // else
                // {
                //     throw new ArgumentException();
                // }
            }
        }

        public char AltLoc { get; private set; }

        private String resName;
        public String ResName
        {
            get
            {
                return resName;
            }
            private set
            {
                resName = value;
                // if (value.Length == 3)
                // {
                //     resName = value;
                // }
                // else
                // {
                //     throw new ArgumentException();
                // }
            }
        }

        public string ChainID { get; private set; }

        public int ResSeq { get; private set; }

        private char iCode;
        public char ICode
        {
            get
            {
                return iCode;
            }
            private set
            {
                if ((value >= 'a' && value <= 'z') || (value >= 'A' && value <= 'Z') || value == ' ')
                {
                    iCode = value;
                }
                else
                {
                    throw new ArgumentException();
                }
            }
        }

        public double X { get; private set; }
        public double Y { get; private set; }
        public double Z { get; private set; }

        public double Occupancy { get; private set; }
        public double TempFactor { get; private set; }

        private String element;
        public String Element
        {
            get
            {
                return element;
            }
            private set
            {
                element = value;
                // if (value.Length == 2)
                // {
                //     element = value;
                // }
                // else
                // {
                //     throw new ArgumentException();
                // }
            }
        }

        private String charge;
        public String Charge
        {
            get
            {
                return charge;
            }
            private set
            {
                if (value.Length == 2)
                {
                    charge = value;
                }
                else
                {
                    throw new ArgumentException();
                }
            }
        }

        public bool IsHetAtm { get; private set; }

        public Atom(int serial, String name, char altLoc, String resName,
            string chainID, int resSeq, char iCode, double x, double y, double z,
            double occupancy, double tempFactor, String element, String charge, bool isHetAtm)
        {
            Serial = serial;
            Name = name.Trim();
            AltLoc = altLoc;
            ResName = resName;
            ChainID = chainID;
            ResSeq = resSeq;
            ICode = iCode;
            X = x;
            Y = y;
            Z = z;
            Occupancy = occupancy;
            TempFactor = tempFactor;
            Element = element;
            Charge = charge;
            IsHetAtm = isHetAtm;
        }

        public Atom(String line)
        {
            if (line.Length != 80)
            {
                throw new FormatException("80 characters expected.");
            }
            else if (!line.Substring(0, 6).Equals("ATOM  ") && !line.Substring(0, 6).Equals("HETATM"))
            {
                throw new FormatException("Prefix \"ATOM  \" or \"HETATM\" expected.");
            }
            else
            {
                int parsedInt;
                double parsedDouble;

                if (Int32.TryParse(line.Substring(6, 5), out parsedInt))
                {
                    Serial = parsedInt;
                }
                else
                {
                    throw new FormatException("serial");
                }

                Name = line.Substring(12, 4).Trim();
                AltLoc = line[16];
                ResName = line.Substring(17, 3);
                ChainID = line.Substring(21, 1);

                if (Int32.TryParse(line.Substring(22, 4), out parsedInt))
                {
                    ResSeq = parsedInt;
                }
                else
                {
                    throw new FormatException("resSeq");
                }

                ICode = line[26];

                if (Double.TryParse(line.Substring(30, 8), NumberStyles.Float, new CultureInfo("en-US"), out parsedDouble))
                {
                    X = parsedDouble;
                }
                else
                {
                    throw new FormatException("x");
                }

                if (Double.TryParse(line.Substring(38, 8), NumberStyles.Float, new CultureInfo("en-US"), out parsedDouble))
                {
                    Y = parsedDouble;
                }
                else
                {
                    throw new FormatException("y");
                }

                if (Double.TryParse(line.Substring(46, 8), NumberStyles.Float, new CultureInfo("en-US"), out parsedDouble))
                {
                    Z = parsedDouble;
                }
                else
                {
                    throw new FormatException("z");
                }

                if (Double.TryParse(line.Substring(54, 6), NumberStyles.Float, new CultureInfo("en-US"), out parsedDouble))
                {
                    Occupancy = parsedDouble;
                }
                else
                {
                    throw new FormatException("occupancy");
                }

                if (Double.TryParse(line.Substring(60, 6), NumberStyles.Float, new CultureInfo("en-US"), out parsedDouble))
                {
                    TempFactor = parsedDouble;
                }
                else
                {
                    throw new FormatException("tempFactor");
                }

                Element = line.Substring(76, 2).Trim();
                Charge = line.Substring(78, 2);

                if (line.Substring(0, 6).Equals("HETATM"))
                {
                    IsHetAtm = true;
                }
                else
                {
                    IsHetAtm = false;
                }
            }
        }

        public override string ToString()
        {
			return (IsHetAtm?"HETATM":"ATOM  ")
                + Serial.ToString().PadLeft(5)
                + " "
                + Name
                + AltLoc
                + ResName
                + " "
                + ChainID
                + ResSeq.ToString().PadLeft(4)
                + ICode
                + "   "
                + X.ToString("0.000", new CultureInfo("en-US")).PadLeft(8)
                + Y.ToString("0.000", new CultureInfo("en-US")).PadLeft(8)
                + Z.ToString("0.000", new CultureInfo("en-US")).PadLeft(8)
                + Occupancy.ToString("0.00", new CultureInfo("en-US")).PadLeft(6)
                + TempFactor.ToString("0.00", new CultureInfo("en-US")).PadLeft(6)
                + "          "
                + Element
                + Charge;

        }

        public class Name_ResSeq_ChainID_AltLoc_Comparer : IEqualityComparer<Atom>
        {
            public bool Equals(Atom a1, Atom a2)
            {
                return a1.Name == a2.Name
                    && a1.ResSeq == a2.ResSeq
                    && a1.ChainID == a2.ChainID
                    && a1.AltLoc == a2.AltLoc;
            }

            public int GetHashCode(Atom a)
            {
                return a.Name.GetHashCode()
                    + a.ResSeq.GetHashCode()
                    + a.ChainID.GetHashCode()
                    + a.AltLoc.GetHashCode();
            }
        }
        
        public int CompareTo(Atom other)
        {
            return this.Serial - other.Serial;
        }
        
		public Point Position(){
			return new Point (X, Y, Z);
		}
		
		public const String NAME_N_AMIDE = "N";
		public const String NAME_H_AMIDE = "H";
		public const String NAME_C_ALPHA = "CA";
		public const String NAME_C_CARB = "C";
		public const String NAME_O_CARB = "O";

		public const String ELEMENT_H = "H";

		public const String CHARGE_ZERO = "  ";

		public bool IsNAmide { get { return Element == "N" && Name == NAME_N_AMIDE; } }
		public bool IsHAmide { get { return Element == "H" && Name == NAME_H_AMIDE; } }
		public bool IsCAlpha { get { return Element == "C" && Name == NAME_C_ALPHA; } }
		public bool IsCCarb { get { return Element == "C" && Name == NAME_C_CARB; } }
		public bool IsOCarb { get { return Element == "O" && Name == NAME_O_CARB; } }

		public bool IsHydrogen{ get { return Element == ELEMENT_H; } }
    }
}

using System;
using System.Collections.Generic;
using System.Linq;

namespace protein
{
	public class OptionParseException : Exception {
		public OptionParseException(String message) : base(message) {}
		public OptionParseException(String message, Exception innerException) : base(message, innerException) {}
	}


	public class Options
	{
		public String GlobalHelp { get; set; }
		public List<Option> OptionList { get; private set; }
		public List<Argument> ArgumentList { get; private set; }
		public const String OPTION_INDENT      = "  ";
		public const String OPTION_HELP_INDENT = "          ";
		public Options ()
		{
			this.OptionList = new List<Option> ();
			this.ArgumentList = new List<Argument> ();
			this.AddOption(new Option (new string[]{ "-h", "--help" }, 0, args => PrintHelp())
				.AddHelp("Print this help message and return 0.")
			);
		}
		public void AddOption(Option option){
			this.OptionList.Add (option);
		}
		public void AddArgument(Argument argument){
			this.ArgumentList.Add (argument);
		}
			
		public bool TryParse(IEnumerable<String> args, out List<String> otherArgs) {
			List<String> argList = new List<String> (args);
			otherArgs = new List<String> ();
			while (argList.Count > 0) {
				String name = argList.Pop ();
				if (name.Length > 0 && name [0] == '-') {
					bool used = false;
					foreach (Option option in OptionList) {
						if (option.Names.Contains (name)) {
							if (argList.Count < option.NumArgs) {
								PrintError ("Option {0} requires {1} arguments.", name, option.NumArgs);
								return false;
							}
							List<String> optionArgs = argList.Pop (option.NumArgs);
							try {
								foreach(var constraint_message in option.Constraints) {
									var constraint = constraint_message.Item1;
									var message = constraint_message.Item2;
									if (!constraint(optionArgs)){
										throw new OptionParseException("\"" + optionArgs.EnumerateWithSeparators(" ") + "\", " + message);
									}
								}
								option.Action (optionArgs);
							} catch (OptionParseException e) {
								PrintError ("Option {0}: {1}", name, e.Message);
								return false;
							}
							used = true;
							break;
						}
					}
					if (!used) { 
						PrintError ("Unknown option: {0}.", name);
						return false;						
					}
				} else {
					otherArgs.Add (name);
				}
			}
			return true;
		}

		public void PrintHelp(){
			Console.WriteLine (GlobalHelp);

			Console.WriteLine ("\nUsage:");
			Console.WriteLine (OPTION_INDENT + "dotnet " + System.AppDomain.CurrentDomain.FriendlyName + ".dll [OPTIONS] " + ArgumentList.Select(a => a.Name).EnumerateWithSeparators(" "));
			Console.WriteLine ("\nArguments:");
			foreach (Argument argument in ArgumentList) {
				String line = OPTION_INDENT + argument.Name;
				//Console.WriteLine (line);
				Lib.WriteInColor (ConsoleColor.Cyan, line + "\n");
				foreach (String help in argument.Helps) {
					Console.WriteLine (OPTION_HELP_INDENT + help);
				}
			}			

			Console.WriteLine ("\nOptions:");
			foreach (Option option in OptionList) {
				String line = OPTION_INDENT + option.Names.EnumerateWithCommas () + " " + option.Parameters.EnumerateWithSeparators (" ");
				//Console.WriteLine (line);
				Lib.WriteInColor (ConsoleColor.Cyan, line + "\n");
				foreach (String help in option.Helps) {
					Console.WriteLine (OPTION_HELP_INDENT + help);
				}
			}
			Environment.Exit (0);
		}

		public static void PrintError(String message, params object[] args){
			Lib.WriteError (String.Format(message, args)
				+ "\nRun   dotnet " + System.AppDomain.CurrentDomain.FriendlyName + ".dll --help   for usage information.");
		}
	}


	public class Option {
		public String[] Names { get; private set; }
		public int NumArgs { get; private set; }
		public Action<IEnumerable<String>> Action;
		public List<(Func<List<String>, bool>, String)> Constraints { get; private set; }
		public List<String> Parameters { get; private set; }
		public List<String> Helps { get; private set; }

		public Option(String[] names, int numArgs, Action<IEnumerable<String>> action){
			if (names.Any (name => name [0] != '-'))
				throw new ArgumentException ("Option names must begin with -.");
			Names = names;
			NumArgs = numArgs;
			Action = action;
			Constraints = new List<(Func<List<String>, bool>, String)> ();
			Parameters = new List<String> ();
			Helps = new List<String> ();
		}

		public Option AddConstraint (Func<List<String>,bool> constraint, String errorMessage){
			Constraints.Add ((constraint, errorMessage));
			return this;
		}

		public Option AddParameter (String parameter){
			Parameters.Add (parameter);
			return this;
		}

		public Option AddHelp (String helpMessage){
			Helps.Add (helpMessage);
			return this;
		}

		public static Option SwitchOption(String[] names, Action<bool> storeResult){
			Action<IEnumerable<String>> action = args => {
				storeResult(true);
			};
			return new Option (names, 0, action);
		}
		public static Option IntOption(String[] names, Action<int> storeResult){
			Action<IEnumerable<String>> action = args => {
				try {
					int i = Int32.Parse (args.First ());
					storeResult(i);
				} catch {
					throw new OptionParseException("Could not parse \"" + args.First () + "\" as integer.");
				}
			};
			return new Option (names, 1, action);
		}
		public static Option DoubleOption(String[] names, Action<double> storeResult){
			Action<IEnumerable<String>> action = args => {
				try {
					double i = Double.Parse (args.First ());
					storeResult(i);
				} catch {
					throw new OptionParseException("Could not parse \"" + args.First () + "\" as float.");
				}
			};
			return new Option (names, 1, action);
		}
		public static Option StringOption(String[] names, Action<String> storeResult){
			Action<IEnumerable<String>> action = args => {
				storeResult(args.First ());
			};
			return new Option (names, 1, action);
		}

		public static Option ChoiceOption(String[] names, Action<String> storeResult, IEnumerable<String> choices){
			Action<IEnumerable<String>> action = args => {
				String choice = args.First ();
				if (choices.Contains (choice)){
					storeResult(choice);
				} else {
					throw new OptionParseException("Invalid choice \"" + choice + "\" (valid choices: " + choices.EnumerateWithCommas () + ").");
				}
			};
			return new Option (names, 1, action);
		}
		public static Option DictionaryChoiceOption<T>(String[] names, Action<T> storeResult, Dictionary<T,String> dictionary){
			Action<IEnumerable<String>> action = args => {
				String choice = args.First ();
				try {
					T value = dictionary.First (kv => kv.Value == choice).Key;
					storeResult(value);
				} catch (InvalidOperationException) {
					throw new OptionParseException("Invalid choice \"" + choice + "\" (valid choices: " + dictionary.Values.EnumerateWithCommas () + ").");
				}
			};
			return new Option (names, 1, action);
		}

		public static Option HelpOption(String[] names, Options options){
			Action<IEnumerable<String>> action = args => {
				options.PrintHelp();
			};
			return new Option (names, 0, action);
		}
	}

	public class Argument {
		public String Name { get; private set; }
		public List<String> Helps { get; private set; }

		public Argument (String name) {
			this.Name = name;
			this.Helps = new List<string> ();
		}

		public Argument AddHelp (String helpMessage){
			Helps.Add (helpMessage);
			return this;
		}

	}

}


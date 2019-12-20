using System;
using System.IO;

using protein.Libraries;
using protein.Json;

namespace protein
{
	public class Config
	{
		public String PymolExecutable { get; private set;}
		public String DsspExecutable { get; private set;}
		public String PymolScriptAlign { get; private set;}
		public String PymolScriptSession { get; private set;}

		public Config (String configFileName) {
			String configFile = Path.Combine (AppDomain.CurrentDomain.BaseDirectory, configFileName);
			if (!File.Exists(configFile)){
				Lib.WriteErrorAndExit("Configuration file \"{0}\" not found.", configFile);
			}
			JsonValue configuration = null;
			try {
				configuration = JsonValue.FromFile(configFile);
			} catch {
				Lib.WriteErrorAndExit ("Error while reading configuration file \"{0}\".", configFile);
			}
			try {
				DsspExecutable = configuration["DsspExecutable"].String;
				PymolExecutable = configuration["PymolExecutable"].String;
				PymolScriptAlign = configuration["PymolScriptAlign"].String;
				PymolScriptSession = configuration["PymolScriptSession"].String;
			} catch (JsonKeyNotFoundException e) {
				Lib.WriteErrorAndExit("Missing value of \"{0}\" in configuration file: {1}", e.MissingKey, configFile);
			} catch {
				Lib.WriteErrorAndExit ("Configuration file has incorrect format: {0}", configFile);
			}
		}
			
	}
}


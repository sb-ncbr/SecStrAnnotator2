
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;

namespace SecStrAnnot2.Cif.Raw
{
    /// <summary>
	/// This mmCIF parser is based on https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax ; 
    /// however, not all details are implemented.
	/// </summary>
    internal class CifParser
    {
        // Constants
        private const bool PRINT_TEXT_IN_COLOR = false; //true;
        private const bool PRINT_TOKENS = false; //true;
        private static int MAX_DIGITS_INT = (int)Math.Floor(Math.Log10(int.MaxValue));
        private static int MAX_DIGITS_LONG = (int)Math.Floor(Math.Log10(long.MaxValue));

        private static char[] ordinaryChars = "!%&()*+,-./0123456789:<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ\\^`abcdefghijklmnopqrstuvwxyz{|}~".ToCharArray();
        private static char[] nonBlankChars = ordinaryChars.Concat("\"#$'_;[]").ToArray();
        private static char[] textLeadChars = ordinaryChars.Concat("\"#$'_ \t[]").ToArray();
        private static char[] anyPrintChars = ordinaryChars.Concat("\"#$'_ \t;[]").ToArray();
        private static char[] whiteSpaceChars = " \t\r\n".ToCharArray();

        private static bool[] isOrdinaryChar = Enumerable.Range(0, 256).Select(c => ordinaryChars.Contains<char>((char)c)).ToArray();
        private static bool[] isNonBlankChar = Enumerable.Range(0, 256).Select(c => nonBlankChars.Contains<char>((char)c)).ToArray();
        private static bool[] isTextLeadChar = Enumerable.Range(0, 256).Select(c => textLeadChars.Contains<char>((char)c)).ToArray();
        private static bool[] isAnyPrintChar = Enumerable.Range(0, 256).Select(c => anyPrintChars.Contains<char>((char)c)).ToArray();
        private static bool[] isWhiteSpaceChar = Enumerable.Range(0, 256).Select(c => whiteSpaceChars.Contains<char>((char)c)).ToArray();


        // Enum types
        private enum TokenType { 
            Data, // DATA_ directive (e.g. data_1tqn)
            Save, // SAVE_ directive (e.g. save_nice_frame, save_)
            Loop, // LOOP_ directive (i.e. loop_)
            Tag, // tag (e.g. _atom_site.id)
            Value, // unquoted value (e.g. ALA, 1.521)
            Single, // single-quoted string (e.g. 'something nice')
            Double,  // double-quoted string (e.g. "something else")
            Semicolon  // semicolon-quoted string (e.g. ;something followed by a newline;)
        };
        
        private enum LexicalState { 
            Out, // out of token (in white-space)
            Comment, // in a comment (e.g. # This is a comment)
            Tag, // in a tag (e.g. _atom_site.id)
            Token, // in a regular unquoted token (value or directive, e.g. ALA, 1.521, data_1tqn)
            Single, // in a single-quoted string (e.g. 'something nice')
            Double,  // in a double-quoted string (e.g. "something else")
            Semicolon  // in a semicolon-quoted string (e.g. ;something followed by a newline;)
        };
        
        private enum SyntacticState {
            Out, // out of data_
            Data, // in data_ (possibly in save_)
            Loop // in a loop_, collecting tags
        }


        // Attributes
        public string Text; // full text of the CIF
        private int textLength; // number of characters of the text
        private int[] start; // position of starting character of the i-th token in the text
        private int[] stop; // position of ending character (exlusive) of the i-th token in the text
        private TokenType[] type; // type of the i-th token
        private int[] blocks; // index of the token which is the i-th DATA_ directive
        private int[] saves; // index of the token which is the i-th SAVE_ directive (both opening and closing)
        private int[] loops; // index of the token which is the i-th LOOP_ directive
        private int[] tags; // index of the token which is the i-th tag

        private int[] firstValueForTag;
        private int[] nValuesForTag;
        public int CountValuesForTag(int iTag) => nValuesForTag[iTag];
        private int[] stepForTag;

        public string[] BlockNames { get; private set; }
        public string[] SaveNames { get; private set; }
        public string[] TagNames { get; private set; }


        internal CifParser(string text){
            SetText(text);
            Program.SetTextDone = DateTime.Now;
            LexicalAnalysis_Goto(); // Lexical analysis is ~ 20% faster with Goto version (1.05 s vs. 1.26 s on 3j3q)
            Program.LexicalDone = DateTime.Now;
            ExtractControlTokenNames();
            Program.ExtractNamesDone = DateTime.Now;
            SyntacticAnalysis();
            Program.SyntacticDone = DateTime.Now;
        }

        private void SetText(string text){
            if (Text != null) {
                throw new InvalidOperationException("Attempting to set text into already initialized " + this.GetType());
            }

            Text = text;
            textLength = Text.Length;

            if (textLength == 0 || !isWhiteSpaceChar[Text[textLength-1]]){
                // append a newline to the end of the text to simplify many things
                Text = Text + "\n";
                textLength = Text.Length;
            }
        }

        private void LexicalAnalysis_NoGoto(){
            List<int> startList = new List<int>(textLength / 4); // guessed number of tokens based on the length of text
            List<int> stopList = new List<int>(textLength / 4); 
            List<TokenType> typeList = new List<TokenType>(textLength / 4);
            List<int> blockIndexList = new List<int>();
            List<int> saveIndexList = new List<int>();
            List<int> loopIndexList = new List<int>(100); // guessed typical number of loops
            List<int> tagIndexList = new List<int>(1000); // guessed typical number of tags
            
            LexicalState state = LexicalState.Out;
            
            for (int i = 0; i < textLength; i++){
                char c = Text[i];
                switch (state){
                    case LexicalState.Out:
                        switch (c){
                            case '#': 
                                state = LexicalState.Comment;
                                startList.Add(i);
                                break;
                            case '\'':
                                state = LexicalState.Single;
                                startList.Add(i);
                                break;
                            case '"':
                                state = LexicalState.Double;
                                startList.Add(i);
                                break;
                            case ';':
                                if (i == 0 || Text[i-1] == '\n' || Text[i-1] =='\r'){
                                    state = LexicalState.Semicolon;
                                    startList.Add(i);
                                } else {
                                    state = LexicalState.Token;
                                    startList.Add(i);
                                }
                                break;
                            default:
                                if (isWhiteSpaceChar[c]){
                                    // stay in Out
                                } else if (isOrdinaryChar[c]){
                                    state = LexicalState.Token;
                                    startList.Add(i);
                                } else if (c == '_'){
                                    state = LexicalState.Tag;
                                    startList.Add(i);
                                } else {
                                    throw NewLexicalCifException(i, state);
                                }
                                break;
                        }
                        break;
                    case LexicalState.Token:
                        if (isNonBlankChar[c]){
                            // stay in Token 
                        } else {
                            state = LexicalState.Out;
                            int s = startList[stopList.Count];
                            if (i - s >= 5 // is this a DATA_ token?
                                && (Text[ s ]=='d' || Text[ s ]=='D') 
                                && (Text[s+1]=='a' || Text[s+1]=='A') 
                                && (Text[s+2]=='t' || Text[s+2]=='T') 
                                && (Text[s+3]=='a' || Text[s+3]=='A') 
                                &&  Text[s+4]=='_') {
                                blockIndexList.Add(stopList.Count);
                                typeList.Add(TokenType.Data);
                                stopList.Add(i);
                            } else if (i - s >= 5  // is this a SAVE_ token?
                                && (Text[ s ]=='s' || Text[ s ]=='S') 
                                && (Text[s+1]=='a' || Text[s+1]=='A') 
                                && (Text[s+2]=='v' || Text[s+2]=='V') 
                                && (Text[s+3]=='e' || Text[s+3]=='E') 
                                &&  Text[s+4]=='_') {
                                saveIndexList.Add(stopList.Count);
                                typeList.Add(TokenType.Save);
                                stopList.Add(i);
                            } else if (i - s == 5  // is this a LOOP_ token?
                                && (Text[ s ]=='l' || Text[ s ]=='L') 
                                && (Text[s+1]=='o' || Text[s+1]=='O') 
                                && (Text[s+2]=='o' || Text[s+2]=='O') 
                                && (Text[s+3]=='p' || Text[s+3]=='P') 
                                &&  Text[s+4]=='_') {
                                loopIndexList.Add(stopList.Count);
                                typeList.Add(TokenType.Loop);
                                stopList.Add(i);
                            } else if (i - s >= 5  // is this a STOP_ token?
                                && (Text[ s ]=='s' || Text[ s ]=='S') 
                                && (Text[s+1]=='t' || Text[s+1]=='T') 
                                && (Text[s+2]=='o' || Text[s+2]=='O') 
                                && (Text[s+3]=='p' || Text[s+3]=='P') 
                                &&  Text[s+4]=='_') {
                                startList.RemoveAt(stopList.Count); // remove from tokens
                            } else if (i - s >= 7  // is this a GLOBAL_ token?
                                && (Text[ s ]=='g' || Text[ s ]=='G') 
                                && (Text[s+1]=='l' || Text[s+1]=='L') 
                                && (Text[s+2]=='o' || Text[s+2]=='O') 
                                && (Text[s+3]=='b' || Text[s+3]=='B') 
                                && (Text[s+4]=='a' || Text[s+4]=='A') 
                                && (Text[s+5]=='l' || Text[s+5]=='L') 
                                &&  Text[s+6]=='_') {
                                startList.RemoveAt(stopList.Count); // remove from tokens
                            } else { // this is a normal value token
                                typeList.Add(TokenType.Value);
                                stopList.Add(i);
                            }
                        }
                        break;
                    case LexicalState.Tag:
                        if (isNonBlankChar[c]){
                            // stay in Token 
                        } else {
                            state = LexicalState.Out;
                            tagIndexList.Add(stopList.Count);
                            typeList.Add(TokenType.Tag);
                            stopList.Add(i);
                        }
                        break;
                    case LexicalState.Single:
                        if (c == '\'' && isWhiteSpaceChar[Text[i+1]]){
                            state = LexicalState.Out;
                            typeList.Add(TokenType.Single);
                            stopList.Add(i+1);
                        } else if (isAnyPrintChar[c]){
                            // stay in Single
                        } else {
                            throw NewLexicalCifException(i, state);
                        }
                        break;
                    case LexicalState.Double:
                        if (c == '"' && isWhiteSpaceChar[Text[i+1]]){
                            state = LexicalState.Out;
                            typeList.Add(TokenType.Double);
                            stopList.Add(i+1);
                        } else if (isAnyPrintChar[c]){
                            // stay in Double
                        } else {
                            throw NewLexicalCifException(i, state);
                        }
                        break;
                    case LexicalState.Semicolon:
                        if (c == ';' &&  (Text[i-1] == '\n' || Text[i-1] =='\r')){
                            state = LexicalState.Out;
                            typeList.Add(TokenType.Semicolon);
                            stopList.Add(i+1);
                        } else {
                            // stay in Semicolon
                        }
                        break;
                    case LexicalState.Comment:
                        if (c == '\n' || c == '\r'){
                            state = LexicalState.Out;
                            startList.RemoveAt(stopList.Count); // remove from tokens
                        } else {
                            // stay in Comment
                        }
                        break;
                    default:
                        throw new NotImplementedException();
                }

                // Print out original text in color
                if (PRINT_TEXT_IN_COLOR){
                    ConsoleColor fgColor;
                    ConsoleColor bgColor;
                    switch(state){
                        case LexicalState.Out: 
                            if (c == '\'' || c == '"') { fgColor = ConsoleColor.Blue; bgColor=ConsoleColor.Gray; }
                            else if (c == ';') { fgColor = ConsoleColor.Green; bgColor=ConsoleColor.Gray; }
                            else { fgColor = ConsoleColor.Black; bgColor=ConsoleColor.Yellow; }
                            break;
                        case LexicalState.Comment: 
                            fgColor = ConsoleColor.Magenta; bgColor=ConsoleColor.Gray; break;
                        case LexicalState.Tag: 
                            fgColor = ConsoleColor.Red; bgColor=ConsoleColor.Gray; break;
                        case LexicalState.Token: 
                            fgColor = ConsoleColor.Black; bgColor=ConsoleColor.Gray; break;
                        case LexicalState.Single:
                        case LexicalState.Double:
                            fgColor = ConsoleColor.Blue; bgColor=ConsoleColor.Gray; break;
                        case LexicalState.Semicolon:
                            fgColor = ConsoleColor.Green; bgColor=ConsoleColor.Gray; break;
                        default:
                            fgColor = ConsoleColor.Black; bgColor=ConsoleColor.Gray; break;
                    }
                    Console.BackgroundColor = bgColor;
                    Console.ForegroundColor = fgColor;
                    Console.Write(c);
                }
            }
            if (state != LexicalState.Out){
                throw new CifException("Unexpected end of file");
            }
            
            // Print out extracted tokens in color
            if (PRINT_TOKENS){
                ConsoleColor origFg = Console.ForegroundColor;
                ConsoleColor origBg = Console.BackgroundColor;
                for (int t = 0; t < startList.Count; t++)
                {
                    if (typeList[t] == TokenType.Tag){
                        Console.WriteLine(" ");
                        Console.BackgroundColor = ConsoleColor.Gray;
                        Console.ForegroundColor = ConsoleColor.Red;
                        Console.Write(Text.Substring(startList[t], stopList[t]-startList[t]));
                        Console.BackgroundColor = ConsoleColor.Yellow;
                        Console.ForegroundColor = ConsoleColor.Black;
                        Console.Write(" ");
                    } else if (typeList[t] == TokenType.Single || typeList[t] == TokenType.Double){
                        Console.BackgroundColor = ConsoleColor.Gray;
                        Console.ForegroundColor = ConsoleColor.Blue;
                        Console.Write(Text.Substring(startList[t], stopList[t]-startList[t]));
                        Console.BackgroundColor = ConsoleColor.Yellow;
                        Console.ForegroundColor = ConsoleColor.Black;
                        Console.Write(" ");
                    } else if (typeList[t] == TokenType.Semicolon){
                        Console.WriteLine(" ");
                        Console.BackgroundColor = ConsoleColor.Gray;
                        Console.ForegroundColor = ConsoleColor.Green;
                        Console.Write(Text.Substring(startList[t], stopList[t]-startList[t]));
                        Console.BackgroundColor = ConsoleColor.Yellow;
                        Console.ForegroundColor = ConsoleColor.Black;
                        Console.Write(" ");
                    } else if (typeList[t] == TokenType.Value){
                        if (typeList[t] == TokenType.Data || typeList[t] == TokenType.Save || typeList[t] == TokenType.Loop) {
                            Console.WriteLine(" ");
                            Console.BackgroundColor = ConsoleColor.Cyan;
                            Console.ForegroundColor = ConsoleColor.Black;
                            Console.Write(Text.Substring(startList[t], stopList[t]-startList[t]));
                            Console.BackgroundColor = ConsoleColor.Yellow;
                            Console.ForegroundColor = ConsoleColor.Black;
                            Console.Write(" ");
                        } else {
                            Console.BackgroundColor = ConsoleColor.Gray;
                            Console.ForegroundColor = ConsoleColor.Black;
                            Console.Write(Text.Substring(startList[t], stopList[t]-startList[t]));
                            Console.BackgroundColor = ConsoleColor.Yellow;
                            Console.ForegroundColor = ConsoleColor.Black;
                            Console.Write(" ");
                        }
                    }
                }
                Console.BackgroundColor = origBg;
                Console.ForegroundColor = origFg;
                Console.WriteLine(" ");
            }

            start = startList.ToArray();
            stop = stopList.ToArray();
            type = typeList.ToArray();
            blocks = blockIndexList.ToArray();
            saves = saveIndexList.ToArray();
            loops = loopIndexList.ToArray();
            tags = tagIndexList.ToArray();
        }

        private void LexicalAnalysis_Goto(){
            List<int> startList = new List<int>(textLength / 4); // guessed number of tokens based on the length of text
            List<int> stopList = new List<int>(textLength / 4); 
            List<TokenType> typeList = new List<TokenType>(textLength / 4);
            List<int> blockIndexList = new List<int>();
            List<int> saveIndexList = new List<int>();
            List<int> loopIndexList = new List<int>(100); // guessed typical number of loops
            List<int> tagIndexList = new List<int>(1000); // guessed typical number of tags
            
            int i = -1;
            char c;

            // implicit goto STATE_OUT;
            
            STATE_OUT: 
                if (++i >= textLength) goto STATE_FINISH;
                c = Text[i];
                switch (c){
                    case '#': 
                        startList.Add(i);
                        goto STATE_COMMENT;
                    case '\'':
                        startList.Add(i);
                        goto STATE_SINGLE;
                    case '"':
                        startList.Add(i);
                        goto STATE_DOUBLE;
                    case ';':
                        if (i == 0 || Text[i-1] == '\n' || Text[i-1] =='\r'){
                            startList.Add(i);
                            goto STATE_SEMICOLON;
                        } else {
                            startList.Add(i);
                            goto STATE_TOKEN;
                        }
                    default:
                        if (isWhiteSpaceChar[c]){
                            goto STATE_OUT;
                        } else if (isOrdinaryChar[c]){
                            startList.Add(i);
                            goto STATE_TOKEN;
                        } else if (c == '_'){
                            startList.Add(i);
                            goto STATE_TAG;
                        } else {
                            throw NewLexicalCifException(i, LexicalState.Out);
                        }
                }

            STATE_COMMENT:
                if (++i >= textLength) throw new CifException("Unexpected end of file");
                c = Text[i];
                if (c == '\n' || c == '\r'){
                    startList.RemoveAt(stopList.Count); // remove from tokens
                    goto STATE_OUT;
                } else {
                    goto STATE_COMMENT;
                }

            STATE_TAG:
                if (++i >= textLength) throw new CifException("Unexpected end of file");
                c = Text[i];
                if (isNonBlankChar[c]){
                    goto STATE_TAG;
                } else {
                    tagIndexList.Add(stopList.Count);
                    typeList.Add(TokenType.Tag);
                    stopList.Add(i);
                    goto STATE_OUT;
                }

            STATE_TOKEN:
                if (++i >= textLength) throw new CifException("Unexpected end of file");
                c = Text[i];
                if (isNonBlankChar[c]){
                    goto STATE_TOKEN; 
                } else {
                    int s = startList[stopList.Count];
                    if (i - s >= 5 // is this a DATA_ token?
                        && (Text[ s ]=='d' || Text[ s ]=='D') 
                        && (Text[s+1]=='a' || Text[s+1]=='A') 
                        && (Text[s+2]=='t' || Text[s+2]=='T') 
                        && (Text[s+3]=='a' || Text[s+3]=='A') 
                        &&  Text[s+4]=='_') {
                        blockIndexList.Add(stopList.Count);
                        typeList.Add(TokenType.Data);
                        stopList.Add(i);
                    } else if (i - s >= 5  // is this a SAVE_ token?
                        && (Text[ s ]=='s' || Text[ s ]=='S') 
                        && (Text[s+1]=='a' || Text[s+1]=='A') 
                        && (Text[s+2]=='v' || Text[s+2]=='V') 
                        && (Text[s+3]=='e' || Text[s+3]=='E') 
                        &&  Text[s+4]=='_') {
                        saveIndexList.Add(stopList.Count);
                        typeList.Add(TokenType.Save);
                        stopList.Add(i);
                    } else if (i - s == 5  // is this a LOOP_ token?
                        && (Text[ s ]=='l' || Text[ s ]=='L') 
                        && (Text[s+1]=='o' || Text[s+1]=='O') 
                        && (Text[s+2]=='o' || Text[s+2]=='O') 
                        && (Text[s+3]=='p' || Text[s+3]=='P') 
                        &&  Text[s+4]=='_') {
                        loopIndexList.Add(stopList.Count);
                        typeList.Add(TokenType.Loop);
                        stopList.Add(i);
                    } else if (i - s >= 5  // is this a STOP_ token?
                        && (Text[ s ]=='s' || Text[ s ]=='S') 
                        && (Text[s+1]=='t' || Text[s+1]=='T') 
                        && (Text[s+2]=='o' || Text[s+2]=='O') 
                        && (Text[s+3]=='p' || Text[s+3]=='P') 
                        &&  Text[s+4]=='_') {
                        startList.RemoveAt(stopList.Count); // remove from tokens
                    } else if (i - s >= 7  // is this a GLOBAL_ token?
                        && (Text[ s ]=='g' || Text[ s ]=='G') 
                        && (Text[s+1]=='l' || Text[s+1]=='L') 
                        && (Text[s+2]=='o' || Text[s+2]=='O') 
                        && (Text[s+3]=='b' || Text[s+3]=='B') 
                        && (Text[s+4]=='a' || Text[s+4]=='A') 
                        && (Text[s+5]=='l' || Text[s+5]=='L') 
                        &&  Text[s+6]=='_') {
                        startList.RemoveAt(stopList.Count); // remove from tokens
                    } else { // this is a normal value token
                        typeList.Add(TokenType.Value);
                        stopList.Add(i);
                    }
                    goto STATE_OUT;
                }

            STATE_SINGLE:
                if (++i >= textLength) throw new CifException("Unexpected end of file");
                c = Text[i];
                if (c == '\'' && isWhiteSpaceChar[Text[i+1]]){
                    typeList.Add(TokenType.Single);
                    stopList.Add(i+1);
                    goto STATE_OUT;
                } else if (isAnyPrintChar[c]){
                    goto STATE_SINGLE;
                } else {
                    throw NewLexicalCifException(i, LexicalState.Single);
                }

            STATE_DOUBLE:
                if (++i >= textLength) throw new CifException("Unexpected end of file");
                c = Text[i];
                if (c == '"' && isWhiteSpaceChar[Text[i+1]]){
                    typeList.Add(TokenType.Double);
                    stopList.Add(i+1);
                    goto STATE_OUT;
                } else if (isAnyPrintChar[c]){
                    goto STATE_DOUBLE;
                } else {
                    throw NewLexicalCifException(i, LexicalState.Double);
                }

            STATE_SEMICOLON:
                if (++i >= textLength) throw new CifException("Unexpected end of file");
                c = Text[i];
                if (c == ';' &&  (Text[i-1] == '\n' || Text[i-1] =='\r')){
                    typeList.Add(TokenType.Semicolon);
                    stopList.Add(i+1);
                    goto STATE_OUT;
                } else {
                    goto STATE_SEMICOLON;
                }
            
            STATE_FINISH:
                // Has processed the whole text, continue in normal execution
            
            // Print out extracted tokens in color
            if (PRINT_TOKENS){
                ConsoleColor origFg = Console.ForegroundColor;
                ConsoleColor origBg = Console.BackgroundColor;
                for (int t = 0; t < startList.Count; t++)
                {
                    if (typeList[t] == TokenType.Tag){
                        Console.WriteLine(" ");
                        Console.BackgroundColor = ConsoleColor.Gray;
                        Console.ForegroundColor = ConsoleColor.Red;
                        Console.Write(Text.Substring(startList[t], stopList[t]-startList[t]));
                        Console.BackgroundColor = ConsoleColor.Yellow;
                        Console.ForegroundColor = ConsoleColor.Black;
                        Console.Write(" ");
                    } else if (typeList[t] == TokenType.Single || typeList[t] == TokenType.Double){
                        Console.BackgroundColor = ConsoleColor.Gray;
                        Console.ForegroundColor = ConsoleColor.Blue;
                        Console.Write(Text.Substring(startList[t], stopList[t]-startList[t]));
                        Console.BackgroundColor = ConsoleColor.Yellow;
                        Console.ForegroundColor = ConsoleColor.Black;
                        Console.Write(" ");
                    } else if (typeList[t] == TokenType.Semicolon){
                        Console.WriteLine(" ");
                        Console.BackgroundColor = ConsoleColor.Gray;
                        Console.ForegroundColor = ConsoleColor.Green;
                        Console.Write(Text.Substring(startList[t], stopList[t]-startList[t]));
                        Console.BackgroundColor = ConsoleColor.Yellow;
                        Console.ForegroundColor = ConsoleColor.Black;
                        Console.Write(" ");
                    } else if (typeList[t] == TokenType.Value){
                        if (typeList[t] == TokenType.Data || typeList[t] == TokenType.Save || typeList[t] == TokenType.Loop) {
                            Console.WriteLine(" ");
                            Console.BackgroundColor = ConsoleColor.Cyan;
                            Console.ForegroundColor = ConsoleColor.Black;
                            Console.Write(Text.Substring(startList[t], stopList[t]-startList[t]));
                            Console.BackgroundColor = ConsoleColor.Yellow;
                            Console.ForegroundColor = ConsoleColor.Black;
                            Console.Write(" ");
                        } else {
                            Console.BackgroundColor = ConsoleColor.Gray;
                            Console.ForegroundColor = ConsoleColor.Black;
                            Console.Write(Text.Substring(startList[t], stopList[t]-startList[t]));
                            Console.BackgroundColor = ConsoleColor.Yellow;
                            Console.ForegroundColor = ConsoleColor.Black;
                            Console.Write(" ");
                        }
                    }
                }
                Console.BackgroundColor = origBg;
                Console.ForegroundColor = origFg;
                Console.WriteLine(" ");
            }

            start = startList.ToArray();
            stop = stopList.ToArray();
            type = typeList.ToArray();
            blocks = blockIndexList.ToArray();
            saves = saveIndexList.ToArray();
            loops = loopIndexList.ToArray();
            tags = tagIndexList.ToArray();
        }

        private string GetTokenString(int iToken){
            return Text.Substring(start[iToken], stop[iToken] - start[iToken]);
        }

        private void ExtractControlTokenNames(){
            int nBlocks = blocks.Length;
            int nSaves = saves.Length / 2; // opening and closing directives => divide by 2
            int nLoops = loops.Length;
            int nTags = tags.Length;

            this.BlockNames = new string[nBlocks];
            this.SaveNames = new string[nSaves];
            this.TagNames = new string[nTags];

            // Tag names
            for (int iTag = 0; iTag < nTags; iTag++)
            {
                int iToken = tags[iTag]; // token index
                TagNames[iTag] = GetTokenString(iToken);
            }

            // Save names
            if (saves.Length % 2 !=0){
                throw new CifException("Missing closing save_ directive");
            }
            for (int iSave = 0; iSave < nSaves; iSave++)
            {
                string opening = GetTokenString(2*iSave);
                string closing = GetTokenString(2*iSave+1);
                if (opening.Length > 5){
                    SaveNames[iSave] = opening.Substring(5);
                } else {
                    throw new CifException("Save frame name must not be empty: " + opening);
                }
                if (closing.Length > 5){
                    throw new CifException("Closing save_ directive cannot be followed by a name: " + closing);
                }
            }

            // Block names
            for (int iBlock = 0; iBlock < nBlocks; iBlock++)
            {
                int iToken = blocks[iBlock]; // token index
                string tokenString = GetTokenString(iToken);
                if (tokenString.Length > 5){
                    BlockNames[iBlock] = tokenString.Substring(5);
                } else {
                    throw new CifException("Data block name must not be empty: " + tokenString);
                }
            }
        }

        private void SyntacticAnalysis(){
            int nBlocks = blocks.Length;
            int nSaves = saves.Length / 2; // opening and closing directives => divide by 2
            int nLoops = loops.Length;
            int nTags = tags.Length;

            int[] firstSaveInBlock = new int[nBlocks];
            //int[] firstLoopInBlock = new int[nBlocks];
            int[] firstTagInBlock = new int[nBlocks];
            int[] firstLoopInSave = new int[nSaves];
            int[] firstTagInSave = new int[nSaves];
            int[] firstLoopAfterSave = new int[nSaves];
            int[] firstTagAfterSave = new int[nSaves];
            int[] firstTagInLoop = new int[nLoops];
            int[] firstTagAfterLoop = new int[nLoops];
            firstValueForTag = new int[nTags];
            nValuesForTag = new int[nTags];
            stepForTag = new int[nTags];

            int[] controlTokens = blocks.Union(saves).Union(loops).Union(tags).Append(start.Length).OrderBy(i=>i).ToArray();
            int nControls = controlTokens.Length - 1;

            int iData = 0; // number of seen data_ tokens
            bool insideSave = false;
            int iSave = 0; // number of seen opening save_ tokens
            int iLoop = 0; // number of seen loop_ tokens
            int iTag = 0; // number of seen tag tokens

            SyntacticState state = SyntacticState.Out;


            for (int iControl = 0; iControl < nControls; iControl++)
            {
                int iToken = controlTokens[iControl];
                TokenType foo = type[iToken]; // rename to tokenType
                switch(state){
                    case SyntacticState.Out:
                        if (foo == TokenType.Data){
                            state = SyntacticState.Data;
                            iData++;
                        } else {
                            throw NewSyntacticCifException(iToken, state, "Expected data_ token");
                        }
                        break;
                    case SyntacticState.Data:
                        switch (foo){
                            case TokenType.Data:
                                // stay in Data, move to the next block
                                firstSaveInBlock[iData] = iSave;
                                firstTagInBlock[iData] = iTag;
                                iData++;
                                break;
                            case TokenType.Save:
                                // stay in Data
                                if (!insideSave){
                                    insideSave = true;
                                    firstLoopInSave[iSave] = iLoop;
                                    firstTagInSave[iSave] = iTag;
                                    iSave++;
                                } else {
                                    insideSave = false;
                                    firstLoopAfterSave[iSave-1] = iLoop;
                                    firstTagAfterSave[iSave-1] = iTag;
                                }
                                break;
                            case TokenType.Loop:
                                state = SyntacticState.Loop;
                                firstTagInLoop[iLoop] = iTag;
                                iLoop++;
                                break;
                            case TokenType.Tag:
                                // stay in Data
                                TokenType nextType = type[iToken+1];
                                if (nextType == TokenType.Value || nextType == TokenType.Single || nextType == TokenType.Double || nextType == TokenType.Semicolon){
                                    firstValueForTag[iTag] = tags[iTag] + 1;
                                    nValuesForTag[iTag] = 1;
                                    stepForTag[iTag] = 0;
                                } else {
                                    throw NewSyntacticCifException(iToken, state, "Missing value for a tag");
                                }
                                iTag++;
                                break;
                        }
                        break;
                    case SyntacticState.Loop:
                        if (foo == TokenType.Tag){
                            TokenType nextType = type[iToken+1];
                            if (nextType != TokenType.Tag){
                                int firstTag = firstTagInLoop[iLoop-1];
                                int nTagsInLoop = iTag - firstTag + 1;
                                int firstValue = iToken + 1;
                                int nValuesInLoop = controlTokens[iControl+1] - firstValue;
                                if (nValuesInLoop % nTagsInLoop != 0){
                                    throw NewSyntacticCifException(loops[iLoop-1], SyntacticState.Loop, "Number of values in a loop must be a multiple of the number of tags");
                                }
                                int nRows = nValuesInLoop / nTagsInLoop;
                                for (int i = 0; i < nTagsInLoop; i++) {
                                    firstValueForTag[firstTag+i] = firstValue + i;
                                    nValuesForTag[firstTag+i] = nRows;
                                    stepForTag[firstTag+i] = nTagsInLoop;
                                }
                                firstTagAfterLoop[iLoop-1] = iTag+1;
                                state = SyntacticState.Data;
                            } // else stay in loop
                            iTag++;
                        } else {
                            throw NewSyntacticCifException(iToken, SyntacticState.Loop, "Unexpected token");
                        }
                        break;
                    default:
                        throw new NotImplementedException();
                }
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal int GetValueStartPosition(int iTag, int iValue){
            int iToken = firstValueForTag[iTag] + iValue * stepForTag[iTag];
            return start[iToken]; 
        }
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal (int, int) GetValuePosition(int iTag, int iValue){
            int iToken = firstValueForTag[iTag] + iValue * stepForTag[iTag];
            return (start[iToken], stop[iToken]); 
        }
        
        internal string GetValueAsString(int iTag, int iValue){
            int iToken = firstValueForTag[iTag] + iValue * stepForTag[iTag];
            return Text.Substring(start[iToken], stop[iToken] - start[iToken]);
        }
        internal string[] GetValuesAsStrings(int iTag){
            int count = nValuesForTag[iTag];
            string[] array = new string[count];
            int iToken = firstValueForTag[iTag];
            int step = stepForTag[iTag];
            for (int iValue = 0; iValue < count; iValue++)
            {
                array[iValue] = Text.Substring(start[iToken], stop[iToken] - start[iToken]);
                iToken += step;
            }
            return array;
        }
        internal string[] GetValuesAsStrings(int iTag, int[] indices){
            int count = indices.Length;
            string[] array = new string[count];
            int first = firstValueForTag[iTag];
            int step = stepForTag[iTag];
            for (int i = 0; i < count; i++)
            {
                int iValue = indices[i];
                int iToken = first + iValue * step;
                array[i] = Text.Substring(start[iToken], stop[iToken] - start[iToken]);
            }
            return array;
        }

        internal char[] GetValuesAsChars(int iTag){
            int count = nValuesForTag[iTag];
            char[] array = new char[count];
            int iToken = firstValueForTag[iTag];
            int step = stepForTag[iTag];
            for (int iValue = 0; iValue < count; iValue++)
            {
                if (stop[iToken] != start[iToken] + 1){
                    throw NewValueParsingCifException(iToken, "Cannot parse as character");
                }
                array[iValue] = Text[start[iToken]];
                iToken += step;
            }
            return array;
        }
        internal char[] GetValuesAsChars(int iTag, int[] indices){
            int count = indices.Length;
            char[] array = new char[count];
            int first = firstValueForTag[iTag];
            int step = stepForTag[iTag];
            for (int i = 0; i < count; i++)
            {
                int iValue = indices[i];
                int iToken = first + iValue * step;
                if (stop[iToken] != start[iToken] + 1){
                    throw NewValueParsingCifException(iToken, "Cannot parse as character");
                }
                array[i] = Text[start[iToken]];
            }
            return array;
        }

        // Assumes that a whitespace will eventually be found in each string (does not check for out-of-bound)
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal bool StringsEqualUntilWhitespace(int start1, int start2) {
            while (true) {
                char c1 = Text[start1];
                char c2 = Text[start2];
                // Lib.WriteLineDebug("    StrEq: " + c1 + " " +c2);
                if (c1 != c2) {
                    // mismatch => check if both strings have ended
                    return isWhiteSpaceChar[c1] && isWhiteSpaceChar[c2];
                } else if (isWhiteSpaceChar[c1]){
                    // both strings have ended, followed by the same whitespace character
                    return true;
                }
                // match without ending
                start1++;
                start2++;
            }
        }

        //TODO make faster by replacing cycle variable checking by checking if character is white-space ?
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal int ParseInteger(int iToken){
            bool negative = false;
            int value = 0;
            int j = start[iToken];
            int stp = stop[iToken];
            // check length
            if (stp - j > MAX_DIGITS_INT){
                try {
                    return int.Parse(Text.Substring(j, stp - j));
                } catch (System.OverflowException e) {
                    throw NewValueParsingCifException(iToken, "Cannot parse as integer (" + e.Message + ")");
                }
            }
            // get sign
            if (Text[j] == '-'){
                negative = true;
                j++;
            } else if (Text[j] == '+'){
                j++;
            }
            // get value
            for ( ; j < stp; j++)
            {
                char c = Text[j];
                if (c < '0' || c > '9'){
                    throw NewValueParsingCifException(iToken, "Cannot parse as integer");
                }
                value *= 10;
                value += (c - '0');
            }
            if (negative){
                value *= -1;
            }
            return value;
        }
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal int ParseInteger(int iToken, int defaultValue){
            bool negative = false;
            int value = 0;
            int j = start[iToken];
            int stp = stop[iToken];
            // check if is not . or ?
            if (j + 1 == stp && (Text[j] == '.' || Text[j] == '?')){
                return defaultValue;
            }
            // check length
            if (stp - j > MAX_DIGITS_INT){
                try {
                    return int.Parse(Text.Substring(j, stp - j));
                } catch (System.OverflowException e) {
                    throw NewValueParsingCifException(iToken, "Cannot parse as integer (" + e.Message + ")");
                }
            }
            // get sign
            if (Text[j] == '-'){
                negative = true;
                j++;
            } else if (Text[j] == '+'){
                j++;
            }
            // get value
            for ( ; j <stp; j++)
            {
                char c = Text[j];
                if (c < '0' || c > '9'){
                    throw NewValueParsingCifException(iToken, "Cannot parse as integer");
                }
                value *= 10;
                value += (c - '0');
            }
            if (negative){
                value *= -1;
            }
            return value;
        }
        internal int[] GetValuesAsIntegers(int iTag){
            int count = nValuesForTag[iTag];
            int[] array = new int[count];
            int iToken = firstValueForTag[iTag];
            int step = stepForTag[iTag];
            for (int iValue = 0; iValue < count; iValue++)
            {
                array[iValue] = ParseInteger(iToken); // parsing manually instead of Int32.Parse() for efficiency
                iToken += step;
            }
            return array;
        }
        internal int[] GetValuesAsIntegers(int iTag, int[] indices){
            int count = indices.Length;
            int[] array = new int[count];
            int first = firstValueForTag[iTag];
            int step = stepForTag[iTag];
            for (int i = 0; i < count; i++)
            {
                int iValue = indices[i];
                int iToken = first + iValue * step;
                array[i] = ParseInteger(iToken); // parsing manually instead of Int32.Parse() for efficiency
            }
            return array;
        }
        internal int[] GetValuesAsIntegers(int iTag, int defaultValue){
            int count = nValuesForTag[iTag];
            int[] array = new int[count];
            int iToken = firstValueForTag[iTag];
            int step = stepForTag[iTag];
            for (int iValue = 0; iValue < count; iValue++)
            {
                array[iValue] = ParseInteger(iToken, defaultValue); // parsing manually instead of Int32.Parse() for efficiency
                iToken += step;
            }
            return array;
        }
        internal int[] GetValuesAsIntegers(int iTag, int[] indices, int defaultValue){
            int count = indices.Length;
            int[] array = new int[count];
            int first = firstValueForTag[iTag];
            int step = stepForTag[iTag];
            for (int i = 0; i < count; i++)
            {
                int iValue = indices[i];
                int iToken = first + iValue * step;
                array[i] = ParseInteger(iToken, defaultValue); // parsing manually instead of Int32.Parse() for efficiency
            }
            return array;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)] // does not work in Debug, but in Release is even faster than manual inlining
        private double ParseDouble(int iToken){
            bool negative = false;
            long whole = 0;
            bool hasFractional = false; // no decimal point
            long fractional = 0;
            long fractionalScale = 1;
            bool hasExponent = false;
            bool negativeExponent = false;
            int exponent = 0;
            double value;
            int j = start[iToken];
            int stp = stop[iToken];
            // check length
            if (stp - j > MAX_DIGITS_LONG){
                try {
                    return double.Parse(Text.Substring(j, stp - j));
                } catch (System.OverflowException e) {
                    throw NewValueParsingCifException(iToken, "Cannot parse as real number (" + e.Message + ")");
                }
            }
            // parse sign
            if (Text[j] == '-'){
                negative = true;
                j++;
            } else if (Text[j] == '+'){
                j++;
            }
            // parse the whole part of the real number             
            for ( ; j < stp; j++)
            {
                char c = Text[j];
                if (c == '.'){
                    hasFractional = true;
                    j++;
                    break;
                }
                if (c < '0' || c > '9'){
                    throw NewValueParsingCifException(iToken, "Cannot parse as real number");
                }
                whole *= 10;
                whole += (c - '0');
            }
            // parse the fractional part of the real number (if any)
            if (hasFractional){
                for ( ; j < stp; j++)
                {
                    char c = Text[j];
                    if (c == 'e' || c == 'E'){
                        hasExponent = true;
                        j++;
                        break;
                    }
                    if (c < '0' || c > '9'){
                        throw NewValueParsingCifException(iToken, "Cannot parse as real number");
                    }
                    fractionalScale *= 10;
                    fractional *= 10;
                    fractional += (c - '0');
                }
                value = whole + (double)fractional / fractionalScale;
            } else {
                value = whole;
            }
            if (negative){
                value *= -1;
            }
            if (hasExponent){
                if (Text[j] == '-'){
                    negativeExponent = true;
                    j++;
                } else if (Text[j] == '+'){
                    j++;
                }
                for ( ; j < stp; j++)
                {
                    char c = Text[j];
                    if (c < '0' || c > '9'){
                        throw NewValueParsingCifException(iToken, "Cannot parse as real number");
                    }
                    exponent *= 10;
                    exponent += (c - '0');
                }
                if (negativeExponent){
                    exponent *= -1;
                }
                value *= Math.Pow(10, exponent);
            }
            return value;
        }
        [MethodImpl(MethodImplOptions.AggressiveInlining)] // does not work in Debug, but in Release is even faster than manual inlining
        private double ParseDouble(int iToken, double defaultValue){
            bool negative = false;
            long whole = 0;
            bool hasFractional = false; // no decimal point
            long fractional = 0;
            long fractionalScale = 1;
            bool hasExponent = false;
            bool negativeExponent = false;
            int exponent = 0;
            double value;
            int j = start[iToken];
            int stp = stop[iToken];
            // check if is not . or ?
            if (j + 1 == stp && (Text[j] == '.' || Text[j] == '?')){
                return defaultValue;
            }
            // check length
            if (stp - j > MAX_DIGITS_LONG){
                try {
                    return double.Parse(Text.Substring(j, stp - j));
                } catch (System.OverflowException e) {
                    throw NewValueParsingCifException(iToken, "Cannot parse as real number (" + e.Message + ")");
                }
            }
            // parse sign
            if (Text[j] == '-'){
                negative = true;
                j++;
            } else if (Text[j] == '+'){
                j++;
            }
            // parse the whole part of the real number             
            for ( ; j < stp; j++)
            {
                char c = Text[j];
                if (c == '.'){
                    hasFractional = true;
                    j++;
                    break;
                }
                if (c < '0' || c > '9'){
                    throw NewValueParsingCifException(iToken, "Cannot parse as real number");
                }
                whole *= 10;
                whole += (c - '0');
            }
            // parse the fractional part of the real number (if any)
            if (hasFractional){
                for ( ; j < stp; j++)
                {
                    char c = Text[j];
                    if (c == 'e' || c == 'E'){
                        hasExponent = true;
                        j++;
                        break;
                    }
                    if (c < '0' || c > '9'){
                        throw NewValueParsingCifException(iToken, "Cannot parse as real number");
                    }
                    fractionalScale *= 10;
                    fractional *= 10;
                    fractional += (c - '0');
                }
                value = whole + (double)fractional / fractionalScale;
            } else {
                value = whole;
            }
            if (negative){
                value *= -1;
            }
            if (hasExponent){
                if (Text[j] == '-'){
                    negativeExponent = true;
                    j++;
                } else if (Text[j] == '+'){
                    j++;
                }
                for ( ; j < stp; j++)
                {
                    char c = Text[j];
                    if (c < '0' || c > '9'){
                        throw NewValueParsingCifException(iToken, "Cannot parse as real number");
                    }
                    exponent *= 10;
                    exponent += (c - '0');
                }
                if (negativeExponent){
                    exponent *= -1;
                }
                value *= Math.Pow(10, exponent);
            }
            return value;
        }
        internal double[] GetValuesAsDoubles(int iTag){
            int count = nValuesForTag[iTag];
            double[] array = new double[count];
            int iToken = firstValueForTag[iTag];
            int step = stepForTag[iTag];
            for (int iValue = 0; iValue < count; iValue++)
            {
                array[iValue] = ParseDouble(iToken); // parsing manually instead of Double.Parse() for efficiency
                iToken += step;
            }
            return array;
        }
        internal double[] GetValuesAsDoubles(int iTag, int[] indices){
            int count = indices.Length;
            double[] array = new double[count];
            int first = firstValueForTag[iTag];
            int step = stepForTag[iTag];
            for (int i = 0; i < count; i++)
            {
                int iValue = indices[i];
                int iToken = first + iValue * step;
                array[i] = ParseDouble(iToken); // parsing manually instead of Double.Parse() for efficiency
            }
            return array;
        }
        internal double[] GetValuesAsDoubles(int iTag, double defaultValue){
            int count = nValuesForTag[iTag];
            double[] array = new double[count];
            int iToken = firstValueForTag[iTag];
            int step = stepForTag[iTag];
            for (int iValue = 0; iValue < count; iValue++)
            {
                array[iValue] = ParseDouble(iToken, defaultValue); // parsing manually instead of Double.Parse() for efficiency
                iToken += step;
            }
            return array;
        }
        internal double[] GetValuesAsDoubles(int iTag, int[] indices, double defaultValue){
            int count = indices.Length;
            double[] array = new double[count];
            int first = firstValueForTag[iTag];
            int step = stepForTag[iTag];
            for (int i = 0; i < count; i++)
            {
                int iValue = indices[i];
                int iToken = first + iValue * step;
                array[i] = ParseDouble(iToken, defaultValue); // parsing manually instead of Double.Parse() for efficiency
            }
            return array;
        }

        internal int[] GetIndicesWhere(int iTag, Func<string,bool> predicate){
            int count = nValuesForTag[iTag];
            List<int> indices = new List<int>();
            int iToken = firstValueForTag[iTag];
            int step = stepForTag[iTag];
            for (int iValue = 0; iValue < count; iValue++)
            {
                string str = Text.Substring(start[iToken], stop[iToken] - start[iToken]);
                if (predicate(str)){
                    indices.Add(iValue);
                }
                iToken += step;
            }
            return indices.ToArray();
        }
        internal int[] GetIndicesWhere(int iTag, Func<string,int,int,bool> predicate){
            int count = nValuesForTag[iTag];
            List<int> indices = new List<int>();
            int iToken = firstValueForTag[iTag];
            int step = stepForTag[iTag];
            for (int iValue = 0; iValue < count; iValue++)
            {
                if (predicate(Text, start[iToken], stop[iToken])){
                    indices.Add(iValue);
                }
                iToken += step;
            }
            return indices.ToArray();
        }

        internal int[] GetIndicesWith(int iTag, string sample){
            int count = nValuesForTag[iTag];
            //string[] array = new string[count];
            int length = sample.Length;
            List<int> indices = new List<int>();
            int iToken = firstValueForTag[iTag];
            int step = stepForTag[iTag];
            for (int iValue = 0; iValue < count; iValue++)
            {
                // checking string equality manually instead of (Text.Substring(...) == sample) for efficiency
                int strt = start[iToken];
                int stp = stop[iToken];
                bool ok = true;
                if (stp-strt == length){
                    for (int i = 0; i < length; i++)
                    {
                        if (Text[strt+i] != sample[i]){
                            ok = false;
                            break;
                        }
                    }
                    if (ok){
                        indices.Add(iValue);
                    }
                }
                iToken += step;
            }
            return indices.ToArray();
        }
        internal int[] GetIndicesWith(int iTag, string[] samples){
            int count = nValuesForTag[iTag];
            //string[] array = new string[count];
            int nSamples = samples.Length;
            int[] lengths = samples.Select(s => s.Length).ToArray();
            List<int> indices = new List<int>();
            int iToken = firstValueForTag[iTag];
            int step = stepForTag[iTag];
            for (int iValue = 0; iValue < count; iValue++)
            {
                // checking string equality manually instead of (Text.Substring(...) == sample) for efficiency
                int strt = start[iToken];
                int stp = stop[iToken];
                for (int s = 0; s < nSamples; s++) // iterate samples
                {
                    string sample = samples[s];
                    int length = lengths[s];
                    bool ok = true;
                    if (stp-strt == length){
                        for (int i = 0; i < length; i++)
                        {
                            if (Text[strt+i] != sample[i]){
                                ok = false;
                                break; // mismatch -> discard sample, continue with next sample
                            }
                        }
                        if (ok){
                            indices.Add(iValue);
                            break; // full match -> discard other samples
                        }
                    }
                }
                iToken += step;
            }
            return indices.ToArray();
        }

        // Gets the indices of values in iValues[from:to], where the value of iTag changes, + to at the end
        // Doesn't check if iValues is non-empty
        //TODO implement variant with int iValueFrom, int iValueTo instead of int[] iValues (faster for sequential iValues)
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private void GetRuns(int iTag, int[] iValues, int from, int to, List<int> outRunStartIndices, List<int> outRunHeadPositions){
            Lib.WriteLineDebug("GetRuns {0} to {1}", from, to);
            int valueCount = to - from;
            int first = firstValueForTag[iTag];
            int step = stepForTag[iTag];
            int iToken = first + iValues[from] * step;
            outRunStartIndices.Clear();
            outRunStartIndices.Add(from);
            int currentHeadPosition = start[iToken]; // position (in Text) of the first token of this run (contiguous region)
            outRunHeadPositions.Clear();
            outRunHeadPositions.Add(currentHeadPosition);
            for (int i = from+1; i < to; i++){
                iToken = first + iValues[i] * step;
                // Lib.WriteLineDebug(GetTokenString(iToken) + " " + StringsEqualUntilWhitespace(start[iToken], currentHeadPosition));
                if (!StringsEqualUntilWhitespace(start[iToken], currentHeadPosition)){
                    // start a new run
                    outRunStartIndices.Add(i);
                    currentHeadPosition = start[iToken];
                    outRunHeadPositions.Add(currentHeadPosition);
                }
            }
            outRunStartIndices.Add(to);
            Lib.WriteLineDebug("run starts: " + outRunStartIndices.Enumerate(" "));
        }
        /// <param name="seekedValue"> position in Text, where the seeked value starts</param>
        /// <param name="values"> positions in Text, where the searched values start</param>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private int LookupValue(List<int> values, int seekedValue, int expectedIndex){
            int nValues = values.Count;
            for (int index = expectedIndex; index < nValues; index++){
                if (StringsEqualUntilWhitespace(seekedValue, values[index])){
                    return index;
                }
            }
            for (int index = 0 ; index < expectedIndex; index++){
                if (StringsEqualUntilWhitespace(seekedValue, values[index])){
                    return index;
                }
            }
            return -1;
        }
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private void GatherSameRuns(int[] array, List<int> startIndicesOfRuns, List<int> headPositionsOfRuns, int[] outArray){
            int nRuns = startIndicesOfRuns.Count - 1; // the last value is sentinel
            List<int> uniqueValueHeadPositions = new List<int>(nRuns); //TODO recycle
            List<int> firstRunForUniqueValues = new List<int>(nRuns);
            List<int> lastRunForUniqueValues = new List<int>(nRuns);
            int[] nextSameRun = new int[nRuns]; // used to implement multiple linked lists  
            int lastUsedUniqueValueIndex = 0;
            for (int iRun = 0; iRun < nRuns; iRun++)
            {
                int uniqueValueIndex = LookupValue(uniqueValueHeadPositions, headPositionsOfRuns[iRun], lastUsedUniqueValueIndex);
                if (uniqueValueIndex == -1){
                    // unique value not found => add a new one
                    uniqueValueHeadPositions.Add(headPositionsOfRuns[iRun]);
                    firstRunForUniqueValues.Add(iRun);
                    lastRunForUniqueValues.Add(iRun);
                } else {
                    // unique value found => add to the correct linked list
                    nextSameRun[lastRunForUniqueValues[uniqueValueIndex]] = iRun;
                    lastRunForUniqueValues[uniqueValueIndex] = iRun; 
                }
            }
            Lib.WriteLineDebug("Unique: " + firstRunForUniqueValues.Enumerate(" "));
            Lib.WriteLineDebug("Linked lists: " + nextSameRun.Enumerate(" "));
            int outCount = startIndicesOfRuns[0];
            for (int iUnique = 0; iUnique < uniqueValueHeadPositions.Count; iUnique++){
                Lib.WriteLineDebug("  unique " + iUnique);
                int iRun = firstRunForUniqueValues[iUnique];
                while (true){
                //for (int iRun = firstRunForUniqueValues[iUnique]; iRun != 0; iRun = nextSameRun[iRun]){ // iterate through the linked list
                    int runStart = startIndicesOfRuns[iRun];
                    int runEnd = startIndicesOfRuns[iRun+1];
                    int runLength = runEnd - runStart;
                    Array.Copy(array, runStart, outArray, outCount, runLength);
                    outCount += runLength;
                    Lib.WriteLineDebug("    run " + iRun);
                    Lib.WriteLineDebug("    run length " + runLength);
                    Lib.WriteLineDebug("    outCount " + outCount);
                    iRun = nextSameRun[iRun];
                    if (iRun == 0){
                        break; // end of linked list
                    }
                }
            }
        }

        internal int[] GatherByValuesInEachRegion(int itag, int[] iValues, int[] startIndicesOfRegions, out int[] startsOfGatheredRuns){
            int nRegions = startIndicesOfRegions.Length - 1; // the last value is sentinel
            int[] gatheredIValues = new int[iValues.Length];
            List<int> recyclableRunStartIndices = new List<int>();
            List<int> recyclableRunHeadPositions = new List<int>();
            for (int iRegion = 0; iRegion < nRegions; iRegion++)
            {
                //TODO put smth. here - call GetRuns and JoinSameRuns (+ change them so that they can with only a part of the iValues array)
                //maybe inline
                GetRuns(itag, iValues, startIndicesOfRegions[iRegion], startIndicesOfRegions[iRegion+1], recyclableRunStartIndices, recyclableRunHeadPositions);
                GatherSameRuns(iValues, recyclableRunStartIndices, recyclableRunHeadPositions, gatheredIValues); //TODO add indexation!!!
                Lib.WriteLineDebug("ivalues: " + iValues.Enumerate(" "));
                Lib.WriteLineDebug("gathered: " + gatheredIValues.Enumerate(" "));
            }
            //TODO fill startsOfGatheredRuns
            startsOfGatheredRuns = null;
            return gatheredIValues;
        }
        
        
        
        private CifException NewLexicalCifException(int position, LexicalState state){
            (int line, int column) = GetLineColumn(position);
            string character;
            if (Text[position] == '\n') 
                character = "line-feed";
            else if (Text[position] == '\r') 
                character = "carriage-return";
            else 
                character = "'" + Text[position] + "'";
            string lineText = new string(Text.Skip(position - (column-1)).TakeWhile(c => c!='\n' && c!='\r').ToArray());
            return new CifException("Lexical error on line " + line + ", column " + column 
                + ", character " + character + " (lexical analyser state: " + state +"):\n" + lineText);
        }
        private CifException NewSyntacticCifException(int token, SyntacticState state, string message){
            (int line, int column) = GetLineColumn(start[token]);            
            return new CifException("Syntactic error on line " + line + ", column " + column 
                + ", token '" + GetTokenString(token) + "' (syntactic analyser state: " + state +"):\n" + message);
        }
        private CifException NewValueParsingCifException(int token, string message){
            (int line, int column) = GetLineColumn(start[token]);            
            return new CifException("Value parsing error on line " + line + ", column " + column 
                + ", token '" + GetTokenString(token) + "':\n" + message);
        }

        private ValueTuple<int,int> GetLineColumn(int position){
            int i = 0;
            int line = 1;
            int column = 1;
            while(i < position){
                char c = Text[i];
                if (c == '\n' || c == '\r'){
                    // New line
                    char next = Text[i+1];
                    if (c == '\r' && next == '\n' || c == '\n' && next == '\r'){
                        // CR-LF or LF-CR => skip 2 characters and start new line
                        line++;
                        column = 1;
                        i += 2;
                    } else {
                        // LF or CR => skip 1 character and start new line
                        line++;
                        column = 1;
                        i++;
                    }
                } else {
                    // Continue line
                    column++;
                    i++;
                }
                
            }
            return new ValueTuple<int,int>(line, column);
        }


    }
}
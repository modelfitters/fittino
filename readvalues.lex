%{
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <y.tab.h>
void yyerror(char*);
extern int yyParseError;
extern int yyNaN;
extern int yyInfinity;
extern int yyInputFileLineNo;
%}

%s newline

WORD             [a-zA-Z][a-zA-Z0-9"+""-""~"".""_"]*
DIGIT            [0-9]
NUMBER           [-+]?{DIGIT}+"."?{DIGIT}*([eE][-+]?{DIGIT}+)?
ENERGYUNIT       [kMGTPE]?eV
CROSSSECTIONUNIT [npfa]b
COMMENT          #.*\n
COMPARATOR       [<>]
GOESTO           ->
BRA              "("
KET              ")"
COMMA            ","
WHITESPACE       [ \t]+
PATH             ["./""/""~"][a-zA-Z0-9"+""-""~"".""_""/"]*
ERRORSIGN        ("+-"|"-+")

%%
{COMMENT}                 {
                              BEGIN(newline);
/*                              printf("Read comment: %s\n", yytext); */
			      sprintf(yylval.name, "%s\n", yytext);
                              return T_NEWLINE;
                          }

BLOCK|Block               {
                              BEGIN(INITIAL);
                              return T_BLOCK;
                          }

DECAY                     {
                              BEGIN(INITIAL);
                              return T_DECAY;
                          }

XS                        {
                              BEGIN(INITIAL);
                              return T_XS;
                          }


Q[ \t]*=                  {
                              BEGIN(INITIAL);
                              return T_SCALE;
                          }

{CROSSSECTIONUNIT}        {
                              if (!strcmp(yytext, "ab")) {
                                  yylval.integer = -3;
                                  return T_CROSSSECTIONUNIT;
                              }
                              else if (!strcmp(yytext, "fb")) {
                                  yylval.integer = 0;
                                  return T_CROSSSECTIONUNIT;
                              }
                              else if (!strcmp(yytext, "pb")) {
                                  yylval.integer = 3;
                                  return T_CROSSSECTIONUNIT;
                              }
                              else if (!strcmp(yytext, "nb")) {
                                  yylval.integer = 6;
                                  return T_CROSSSECTIONUNIT;
                              }
                              else {
			        fprintf(stderr, "Unknown cross section unit in input file line %d", yyInputFileLineNo);
				yyParseError = 1;
				/* yyerror("Unknown cross section unit"); */
                              }
                          }
 
{GOESTO}                  {
                              return T_GOESTO;
                          }
 
 
{COMPARATOR}              {
                              strcpy(yylval.name, yytext);
                              return T_COMPARATOR;
                          }
 
{COMMA}                   {
                              return T_COMMA;
                          }
 
{BRA}                     {
                              return T_BRA;
                          }

{KET}                     {
                              return T_KET;
                          }
 
<newline>{WORD}           {
                              BEGIN(INITIAL);
                              if (strlen(yytext) > 254) {
				fprintf(stderr, "Key too long in input file line %d", yyInputFileLineNo);
				yyParseError = 1;
                                /*  yyerror("Key too long"); */
                              }
			      if (!strcmp(yytext, "universality")) { 
				  return T_UNIVERSALITY;
                              }
			      if (!strcmp(yytext, "Calculator")) { 
				  return T_CALCULATOR;
                              }
			      if (!strcmp(yytext, "RelicDensityCalculator")) { 
				  return T_RELICDENSITYCALCULATOR;
			      }
                              if (!strcmp(yytext, "LEOCalculator")) { 
				  return T_LEOCALCULATOR;
			      }
                              if (!strcmp(yytext, "MarkovInterfaceFilePath")) { 
				  return T_MARKOVINTERFACEFILEPATH;
			      }
			      if (!strcmp(yytext, "BR")) {
				  return T_BR;
			      }
			      if (!strcmp(yytext, "LEObs")) {
				  return T_LEO;
			      }
			      if (!strcmp(yytext, "nofit")) {
				  return T_NOFIT;
			      }
			      if (!strcmp(yytext, "nofitLEObs")) {
				  return T_NOFITLEO;
			      }
                              strcpy(yylval.name, yytext);
                              return T_KEY;
                          }

{ENERGYUNIT}              {
                              if (!strcmp(yytext, "eV")) {
                                  yylval.integer = -9;
                                  return T_ENERGYUNIT;
                              }
                              else if (!strcmp(yytext, "keV")) {
                                  yylval.integer = -6;
                                  return T_ENERGYUNIT;
                              }
                              else if (!strcmp(yytext, "MeV")) {
                                  yylval.integer = -3;
                                  return T_ENERGYUNIT;
                              }
                              else if (!strcmp(yytext, "GeV")) {
                                  yylval.integer = 0;
                                  return T_ENERGYUNIT;
                              }
                              else if (!strcmp(yytext, "TeV")) {
                                  yylval.integer = 3;
                                  return T_ENERGYUNIT;
                              }
                              else if (!strcmp(yytext, "PeV")) {
                                  yylval.integer = 6;
                                  return T_ENERGYUNIT;
                              }
                              else if (!strcmp(yytext, "EeV")) {
                                  yylval.integer = 9;
                                  return T_ENERGYUNIT;
                              }
                              else {
			        fprintf(stderr, "Unknown energy unit in input file line %d", yyInputFileLineNo);
				yyParseError = 1;
                                /*  yyerror("Unknown energy unit"); */
                              }
                          }
 
{ERRORSIGN}               {
                              return T_ERRORSIGN;
                          }

{NUMBER}|"++++++++++++++++"|"????????????????"|"NaN"|"+Inf"|"-Inf" { /* Number, SPheno junk, NaN or infinity */
                              BEGIN(INITIAL);
                              if (!strcmp(yytext, "++++++++++++++++") ||
				  !strcmp(yytext, "????????????????")) {
				  yyParseError = 1;
                                  yylval.real = 0;
			      }
			      else if (!strcmp(yytext, "NaN")) {
				  yyNaN = 1;
                                  yylval.real = 0;
			      }
			      else if (!strcmp(yytext, "+Inf") ||
				       !strcmp(yytext, "-Inf")) {
				  yyInfinity = 1;
                                  yylval.real = 0;
			      }
                              else yylval.real = atof(yytext);
                              return T_NUMBER;
                          }

on|off                    {
                              if (!strcmp(yytext, "on"))       yylval.integer = 1;
                              else if (!strcmp(yytext, "off")) yylval.integer = 0;
                              else               
				{
				  fprintf(stderr, "Unknown switch state in input file line %d", yyInputFileLineNo);
				  yyParseError = 1;
				  /* yyerror("Unknown switch state"); */ }
                              return T_SWITCHSTATE;
                          }
 
alias                     {
                              return T_ALIAS;
                          }

scaling                   {
                              return T_SCALING;
                          }

{WORD}                    {
                              BEGIN(INITIAL);
                              if (strlen(yytext) > 254) {
				yyParseError = 1;
				fprintf(stderr, "Path too long in input file line %d", yyInputFileLineNo);
	                        /*  yyerror("Word too long"); */
			      }
                              strcpy(yylval.name, yytext);
                              return T_WORD;
                          }

{PATH}                    {
                              BEGIN(INITIAL);
                              if (strlen(yytext) > 254) {
				yyParseError = 1;
				fprintf(stderr, "Path too long in input file line %d", yyInputFileLineNo);
	                        /*  yyerror("Path too long"); */
			      }
                              strcpy(yylval.name, yytext);
                              return T_PATH;
                          }

\n                        {
                              BEGIN(newline);
			      strcpy(yylval.name, "\n");
                              return T_NEWLINE;
                          }

{WHITESPACE}              /* Eat up whitespace */

.                         {
                              fprintf(stderr, "Unrecognized character in input file line %d: %s\n", yyInputFileLineNo, yytext);
			      yyParseError = 1;
                              /* yyerror("Unrecognized character"); */
                          }

%%
 
int yywrap(void) {
    return 1;
}

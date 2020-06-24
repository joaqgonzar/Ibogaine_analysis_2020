%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AUTHOR:         Nicolas Rubido                                       %%
%%% CONTACT:        nrubido@fisica.edu.uy | s01nr0@abdn.ac.uk            %%
%%% AFFILIATIONS:                                                        %%
%%%     1 - Universidad de la Republica (UdelaR), Instituto de Fisica    %%
%%%         de Facutad de Ciencias (IFFC), Montevideo 11400, Uruguay.    %%
%%%     2 - University of Aberdeen (UoA), Aberdeen Biomedical Imaging    %%
%%%         Centre (ABIC), Aberdeen AB25 2ZG, United Kingdom.            %%
%%%     3 - University of Aberdeen (UoA), Institute for Complex Systems  %%
%%%         and Mathematical Biology (ICSMB), Aberdeen AB24 3UE, U.K.    %%
%%% VERSION:        24/06/2020                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GENERAL HELP:
%   This function creates an Ordinal Pattern (OP) encoding of a real-valued
% input time-series. The encoding follows the OP definition and methodology
% proposed by C. Bandt and B. Pompe in: Phys. Rev. Lett. 88, 174102 (2002).
% Specifically, the function encodes input signal subsets of length EMBDIM,
% namely, words, into integer-valued symbols, namely, OPs. The resultant OP
% of any word is the ordering needed to set all values in an increasing way
% such that the resultant symbol is the number of permutations needed. This
% encoding produces a symbolic sequence, namely, an OP output series. These
% symbols are then used to find a sliding window information average, i.e.,
% the Shannon OP entropy of the input signal sliding window + its variance.
% As a comparisson, the classic Shannon entropy from bins is also computed.
% Also, the Renyi entropy is computed from the OPs to have further insights
% about the input's complexity and randomness.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ENTROP,PROBSW] = TimeSeries2OPentropy(EEG)
%
global EMBDIM DELAYS


INDATA = EEG; 
SMALLI = 1E-6; 
RAND_T = rand(length(INDATA),1);    % uniformly distributed random series
AUXVEC = SMALLI*(2*RAND_T - 1); % small-intensity white-noise signal
INDATA = INDATA + AUXVEC;       % lift input degeneracies with noise

%
%%% SYMBOLIC ENCODING PARAMETERS:
MINVAL = 5E3;                   % minimum number of symbols per bin or OP
EMBDIM = 3;                     % embedding dimension for OPs
DELAYS = EMBDIM;                % overlapping distance between OPs
NPOSIB = factorial(EMBDIM);     % number of possible OPs for EMBDIM
N_BINS = EMBDIM*NPOSIB;         % number of bins used per sliding window
WINDOW = N_BINS*MINVAL;         % number of data points per sliding window
NUM_QS = 3;			% different q-values for Renyi entropy

%_________________________________________________________________________%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MAIN  PROGRAM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%% DATA INFORMATION ANALYSIS:
NTIMES = cell(1,1);        % time-series length memory allocation
PROBSW = cell(1,2);        % OP + bin probabilities memory allocation
ENTROP = cell(1,3);        % all entropies memory allocation
                           % ...start analysis of sleep states...
    %
    % SLIDING-WINDOW DATA ANALYSIS
    NTIMES = length(INDATA);% sleep-state time-series length
    AUXNUM = NTIMES/WINDOW;     % approx. number of sliding windows
    NUM_SE = floor(AUXNUM);         % minimum number of sliding windows
    BIN_SE = zeros(NUM_SE,2);       % bin Shannon-entropy memory allocation
    OPS_SE = zeros(NUM_SE,2);       % OPs Shannon-entropy memory allocation
    OPS_RE = zeros(NUM_SE,NUM_QS);  % OPs Renyi-entropy memory allocation
    ENDING = NTIMES -(WINDOW-1); % ending index for analysis
    %
    COUNTS = 1;                     % initialise window counter
    BINMAT = zeros(N_BINS,NUM_SE);  % bin prob. matrix memory allocation
    OP_MAT = zeros(NPOSIB,NUM_SE);  % OPs prob. matrix memory allocation
    for TT = 1:WINDOW:ENDING      % ...start sliding-window analysis...
        SLIDEW = TT + WINDOW - 1;       % sliding-window ending index
        VECTOR = INDATA(TT:SLIDEW); % sleep-state window vector
        %
        % HISTOGRAM-BIN ENCODING
        NFREQS = hist( VECTOR, N_BINS );% time-series appearence frequency 
        BINPDF = NFREQS'/sum(NFREQS);   % probability distribution function
        BINMAT(:,COUNTS) = BINPDF;      % allocating bin probabilitiies 
        %
        % SHANNON ENTROPY FOR BINS
        AUXVAL = ShannonEntropy(BINPDF,N_BINS);
        BIN_SE(COUNTS,:) = AUXVAL;      % bin encoding Shannon entropy
        %
        % ORDINAL PATTERN ENCODING
        OPCODE = Signal2OrdPats(VECTOR);% OP encoded time-series
        NFREQS = hist(OPCODE, 1:NPOSIB);% OP appearence frequency
        OPPROB = NFREQS'/sum(NFREQS);   % probability of having the OPs
        OP_MAT(:,COUNTS) = OPPROB;      % allocating OP probabilities
        %
        % SHANNON ENTROPY FOR OPs
        AUXVAL = ShannonEntropy(OPPROB,NPOSIB);
        OPS_SE(COUNTS,:) = AUXVAL;      % allocating OP Shannon entropy
        %
        % RENYI ENTROPY FOR OPs
        AUXVEC = RenyiEntropy(OPCODE,NPOSIB,NUM_QS);
        OPS_RE(COUNTS,:) = AUXVEC;      % allocating OPs Renyi entropies
        %
        COUNTS = COUNTS + 1;            % update window counter
        %
    end;                            % ...end sliding-window analysis...
    %
    PROBSW{1} = BINMAT;          % sliding-window bin probabilities
    PROBSW{2} = OP_MAT;          % sliding-window OP probabilities
    ENTROP{1} = BIN_SE;          % allocation of bin Shannon entropy
    ENTROP{2} = OPS_SE;          % allocation of OPs Shannon entropy
    ENTROP{3} = OPS_RE;          % allocation of OPs Renyi entropy
    %
    disp( [' ... ', '-state analysis (OPs with D = ', ...
          num2str(EMBDIM), ') ready ...']);
    %
                               % ...end analysis of sleep states...
%
%%% SAVE RESULTS/OUTPUT:
return;
end
%%
%
%_________________________________________________________________________%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SUBROTUINES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%_________________________________________________________________________%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function ODATAS = Signal2OrdPats(IDATAS)
%
%  This subroutine encodes a real-valued time-series into an integer-valued
% series, namely, a set of Ordinal Patterns (OPs). The encoding is obtained
% from the subroutine "Permut_Order", which compares the values within each
% sub-set of the time-series that "Signal2OrdPats" provides as input, i.e.,
% words of length EMBDIM. The resultant OP output from "Permut_Order" comes
% from a lexicographic dictionary that is based on ordering the word into a
% increasing-valued word by performing permutations, as explained within.
%
global  EMBDIM  DELAYS
%
N_DATA = length(IDATAS);        % length of input data set
ODATAS = zeros(1,2);            % OP output initialisation
%
%%% TIMES-SERIES ENCODING INTO ORDINAL PATTERNS (OPs)
JJ = 1;                         % initialize allocation counter
ENDING = N_DATA - (EMBDIM - 1); % ending index for the encoding
for NTAU = 1:DELAYS:ENDING,     % ... start encoding loop ...
    INDEXS = NTAU:NTAU+(EMBDIM-1);  % time-series indexes for OP symbol
    T_WORD = IDATAS(INDEXS);        % word formed from input vaules
    PORDER = Permut_Order(T_WORD);  % resultant OP symbol for the word
    ODATAS(JJ) = PORDER;            % OP allocation into output
    JJ = JJ + 1;                    % update allocation counter
end;                            % ... end encoding loop ...
%
return;
end
%
%_________________________________________________________________________%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function PORDER = Permut_Order(T_WORD)
%
%   This subroutine performs a series of permutations to arrange the input,
% T_WORD, into a monotonically increasing word of size: EMBDIM. The output,
% PORDER, is the number of permutations needed to create this ordering + 1.
%   This choice of encoding defines a lexicographic dictionary as a result.
% For example, if EMBDIM = 3 (T_WORD has 3 consecutive time-series values),
% the six possible output values for PORDER are:
%           WORDS	 PORDER
%           (123)       1
%           (132)       2
%           (213)       3
%           (231)       4
%           (312)       5
%           (321)       6
% where, (123), indicates that the difference between consectuvie values is
% increasing, contrary to (321). In other words, (123), requires no further
% ordering, hence, zero permutations (i.e., PORDER = 1). In opposition, the
% word (321) requires 5 permutations in total (i.e., PORDER = 3 + 2 + 1).
%
%%% PERMUTATION ORDER FOR EMBEDING WORD:
AUXVAR = T_WORD;                % auxiliary word allocation
STATUS = 0;                     % initialize the number of permutations
for NN = 1:length(T_WORD)-1,    % ... start word ordering loop ...
    %
    % FIND MINIMUM VALUE OF THE WORD
    MINVAL = min(AUXVAR);           % minimum (variable) word value
	CHOOSE = find(AUXVAR == MINVAL);% location of minimum
    W_LONG = length(AUXVAR);        % length of (variable) word
	POSIBS = factorial(W_LONG - 1); % remaining permutations in word
    %
    % UPDATE THE STATUS OF THE WORD WITHIN THE LEXICOGRAPH ASSIGNED
    UPDATE = (CHOOSE-1)*POSIBS;     % updated state of the permutations
	STATUS = STATUS + UPDATE;       % update symbol of OP (permutations)
    INDEXS = 1:W_LONG;              % indexes of word values
    NEWIND = INDEXS ~= CHOOSE;      % remove minimum from word
    if sum(NEWIND) >= 1,            % ... if there are more letters ...
        TEMPOR = AUXVAR(NEWIND);        % remove used letters in encoding
        AUXVAR = TEMPOR;                % update word without former min.
    else                            % ... if all letters have been used ...
        disp('fatal error!!! verify words conditionals trueness...');
    end;                            % ... end conditionals
end;                            % ... end word ordering loop ...
%
%%% RESULTING PERMUTATION
PORDER = STATUS + 1;            % final permutation status
%
return;
end
%
%_________________________________________________________________________%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [SH_ENT,SH_VAR] = ShannonEntropy(INPROB,N_SYMB)
%
%   This subrotuine computes the Shannon entropy [C. E. Shannon, Bell Syst.
% Tech. J. 27, 623-656 (1948)] for the input probabilities, INPROB, and the
% number of symbols used, N_SYMB, in the alphabet. INPROBs must be a column
% vector. Also, the variance of this quantity is computed, namely, SH_VAR.
%
NONULL = INPROB > 0;            % index for non-null probabilities
NEW_PS = INPROB(NONULL);        % new vector with only non-null probs.
AUXVAL = (NEW_PS')*log2(NEW_PS);% - Shannon entropy in bits per symbol
NEWVAL = -AUXVAL/log2(N_SYMB);	% Shannon entropy in N_SYMB per symbol
SH_ENT = NEWVAL;                % Shannon entropy output allocation
%
AUXVAL = (NEW_PS')*(log2(NEW_PS).^2);
NEWVAL = -AUXVAL/log2(N_SYMB)^2;% quadratic entropy in N_SYMB per symbol
SH_VAR = NEWVAL - SH_ENT^2;     % Shannon entropy variance
%
return;
end
%
%_________________________________________________________________________%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function RENENT = RenyiEntropy(INPUTS,N_SYMB,NUM_QS)
%
%    This subrotuine computes the Rnyi entropy (RENENT) as in [T. Van Erven
% and P. Harremoes, IEEE Trans. Inf. Theor. 60(7), 3797-3820 (2014)] for an
% input series of ordinal patterns, INPUTS, with N_SYMB different patterns.
% For example, q = 2 Renyi-entropy is the knwon Collision entropy (i.e., an
% entropy for independent and identically distributed variables); being the
% smallest integer q-value that this function deals with.
%
NFREQS = hist(INPUTS,1:N_SYMB); % ordinal pattern appearence frequency
OPSPDF = NFREQS/sum(NFREQS);    % OP probability density function
NONULL = OPSPDF > 0;            % index for non-null probabilities
NEW_PS = OPSPDF(NONULL);        % new vector with only non-null probs.
RENENT = zeros(1,NUM_QS);       % Renyi entropy output memory allocation
for QQ = 1:NUM_QS,              % ...start Renyi q-values loop...
    POWERP = NEW_PS.^(QQ+1);        % probability powers
    AUXVAL = log(sum(POWERP))/(-QQ);% Renyi entropy of q-value
    NEWVAL = AUXVAL/log(N_SYMB);	% Renyi entropy in N_SYMB per symbol
    RENENT(QQ) = NEWVAL;            % Renyi entropy output allocation
end;                            % ...end Renyi q-values loop...
%
return;
end
%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%.........................................................................%
%... oooO ............      Programmer footprint:       ~^,        .......%
%.. (   ) .. Oooo ....         N1C07A5 RU81D0          |   \^,     .......%
%... \ ( ... (   ) ...                                 |*     \_   .......%
%.... \_) ... ) / ....           June, 2020           ,^        }  .......%
%........... (_/ .....   FREE SOFTWARE DISTRIBUTION   |        _>  .......%
%....................................................  \_-_ _??/  ........%
%..........................................................*..............%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

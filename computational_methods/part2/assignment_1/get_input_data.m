function [x, Tnod, Tmat, Tdof, Tdn] = get_input_data()

% Node coordinates
x = [%   X       Y       Z
     0.000,  0.000,  0.000; % 1
    -0.725, -0.731,  3.250; % 2
     0.725, -0.731,  3.250; % 3
     0.725,  0.731,  3.250; % 4
    -0.725,  0.731,  3.250; % 5
     2.900, -1.450,  5.630; % 6
     2.900,  0.000,  5.820; % 7
     2.900,  1.450,  5.630; % 8
     0.000,  1.450,  6.340; % 9
     0.000,  0.000,  6.530; % 10
     0.000, -1.450,  6.340; % 11
    -2.900,  1.450,  5.630; % 12
    -2.900,  0.000,  5.820; % 13
    -2.900, -1.450,  5.630; % 14
];

% Bar connectivities
Tnod = [%   A     B
         1     2  % 1
         1     3  % 2
         1     4  % 3
         1     5  % 4
         2     3  % 5
         3     4  % 6
         4     5  % 7
         5     2  % 8
         2     4  % 9
         2    10  % 10
         2    11  % 11
         2    13  % 12
         2    14  % 13
         3     6  % 14
         3     7  % 15
         3    10  % 16
         3    11  % 17
         4     7  % 18
         4     8  % 19
         4     9  % 20
         4    10  % 21
         5     9  % 22
         5    10  % 23
         5    12  % 24
         5    13  % 25
        14    11  % 26
        11    10  % 27
        14    10  % 28
        14    13  % 29
        13    10  % 30
        13    12  % 31
        12    10  % 32
        12     9  % 33
         9    10  % 34
         9     8  % 35
         8     7  % 36
        10     7  % 37
        10     8  % 38
        11     6  % 39
         6     7  % 40
        10     6  % 41
];

% Material connectivities
Tmat = [% Mat. index
         1  % 1
         1  % 2
         1  % 3
         1  % 4
         2  % 5
         2  % 6
         2  % 7
         2  % 8
         2  % 9
         1  % 10
         1  % 11
         1  % 12
         1  % 13
         1  % 14
         1  % 15
         1  % 16
         1  % 17
         1  % 18
         1  % 19
         1  % 20
         1  % 21
         1  % 22
         1  % 23
         1  % 24
         1  % 25
         2  % 26
         2  % 27
         2  % 28
         2  % 29
         2  % 30
         2  % 31
         2  % 32
         2  % 33
         2  % 34
         2  % 35
         2  % 36
         2  % 37
         2  % 38
         2  % 39
         2  % 40
         2  % 41
];

% Degrees of freedom conectivities
Tdof = [%    A            B
         1   2   3   4    5   6  % 1
         1   2   3   7    8   9  % 2
         1   2   3   10   11  12  % 3
         1   2   3   13   14  15  % 4
         4   5   6   7    8   9  % 5
         7   8   9   10   11  12  % 6 
         10  11  12  13   14  15  % 7 
         13  14  15  4    5   6  % 8
         4   5   6   10   11  12  % 9
         4   5   6   28   29  30  % 10
         4   5   6   31   32  33  % 11
         4   5   6   37   38  39  % 12
         4   5   6   40   41  42  % 13
         7   8   9   16   17  18  % 14
         7   8   9   19   20  21  % 15
         7   8   9   28   29  30  % 16
         7   8   9   31   32  33  % 17
         10  11  12  19   20  21  % 18
         10  11  12  22   23  24  % 19
         10  11  12  25   26  27  % 20
         10  11  12  28   29  30  % 21
         13  14  15  25   26  27  % 22
         13  14  15  28   29  30  % 23
         13  14  15  34   35  36  % 24
         13  14  15  37   38  39  % 25
         40  41  42  31   32  33  % 26
         31  32  33  28   29  30  % 27
         40  41  42  28   29  30  % 28
         40  41  42  37   38  39  % 29
         37  38  39  28   29  30  % 30
         37  38  39  34   35  36  % 31
         34  35  36  28   29  30  % 32
         34  35  36  25   26  27  % 33
         25  26  27  28   29  30  % 34
         25  26  27  22   23  24  % 35
         22  23  24  19   20  21  % 36
         28  29  30  19   20  21  % 37
         28  29  30  22   23  24  % 38
         31  32  33  16   17  18  % 39
         16  17  18  19   20  21  % 40
         28  29  30  16   17  18  % 41
];

%DOF per node
Tdn=[1 2 3; %node 1
    4 5 6; %node 2
    7 8 9; %node 3
    10 11 12; %node 4
    13 14 15; %node 5
    16 17 18; %node 6
    19 20 21; %node 7
    22 23 24; %node 8
    25 26 27; %node 9
    28 29 30; %node 10
    31 32 33; %node 11
    34 35 36; %node 12
    37 38 39; %node 13
    40 41 42; %node 14
];


end
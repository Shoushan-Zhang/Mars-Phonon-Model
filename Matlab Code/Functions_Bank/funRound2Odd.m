function oddNumber = funRound2Odd(num)
    roundedNum = round(num);  % Round to the nearest integer
    if mod(roundedNum, 2) == 0
        % If the rounded number is even, add 1 to make it odd
        if num > roundedNum
            oddNumber = roundedNum + 1;
        else
            oddNumber = roundedNum - 1;
        end
    else
        % If already odd, no change needed
        oddNumber = roundedNum;
    end
end

try
    pa = modpar26(1);
    tdm26('tbabr', pa);
catch e
    fprintf('Error: %s\n', e.message);
end
exit;

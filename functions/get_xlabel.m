function xlabel = get_xlabel(j0)

    [jd, jdfrac] = whole_and_frac(j0);

    [year, mon, day, hr, min, sec] = invjday( jd, jdfrac );

    xlabel = ['Days from ', num2str(year), '.', double_digit(mon), '.', ...
              double_digit(day), ' ', double_digit(hr), ':', ...
              double_digit(min), ':', double_digit(sec)];

end

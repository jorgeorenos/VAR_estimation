function Quarterly(df::DataFrame, date::Symbol; drop_incomplete::Bool=true)
    # 1. year and quarter columns
    years      = year.(df[!, date])
    quarter = floor.((month.(df[!, date]) .- 1) ./ 3) .+ 1

    # 2. Create a copy of the dataframe with the new columns
    df2 = copy(df)
    df2[!,:year]       = years
    df2[!,:quarter] = quarter

    # 3. numeric columns
    cols_num = filter(col -> eltype(df2[!, col]) <: Number, names(df2))

    # 4. group by year and quarter, then aggregate by mean
    g = groupby(df2, [:year, :quarter])
    # add a count column to identify incomplete quarters
    agg = combine(g,
                  cols_num .=> mean .=> cols_num,
                  nrow     .=> :count)

    # 5. Drop incomplete quarters if specified
    if drop_incomplete
        agg = filter(row -> row.count == 3, agg)
    end

    # 6. Create a date column for the start of the quarter
    agg[!,:dates] = Date.(agg.year,
                           (agg.quarter .- 1) .* 3 .+ 1,
                           1)

    # 7. Ordenar y seleccionar sólo date+variables
    sort!(agg, :dates)
    select(agg, :dates, cols_num...)
end
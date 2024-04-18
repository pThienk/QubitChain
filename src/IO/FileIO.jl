# Contains functions for abstract File input/output

function interface_IO(filename::String, type::Symbol, msg...)

    if type == :open_write
        write(filename, "")
        write_file(filename, msg)
    elseif type == :open_read

    end
end

function read_file(stream::IOStream)
    # Not used for the moment
end

function write_file(filename::String, msg...)

    file::IOStream = open(filename, "w")

    write(file, join(msg, "\n"))

    close(file)
end
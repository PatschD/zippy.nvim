local M = {}

local ts_utils = require("nvim-treesitter.ts_utils")

local lua_keys = {
  nodes = {
    assignment_statement = true,
    function_call = true,
    arguments = true,
    parameters = true,
    condition = true,
  },
  options = {
    language_format_behaviour = "range",
    print_function = "print",
  },
}

local js_keys = {
  nodes = {
    variable_declarator = true,
    formal_parameters = "statement_block",
    parenthesized_expression = true,
    expression_statement = true,
    identifier = "statement_block",
  },
  options = {
    sibling_block_placement_behaviour = "inside",
    language_format_behaviour = "range",
    print_function = "console.log",
    format_fn = function(t)
      return '"'
          .. "file: "
          .. t.filename
          .. "~"
          .. "line: "
          .. t.line_nr
          .. t.breadcrumbs
          .. "~"
          .. t.current_text
          .. '", '
          .. t.current_text
    end,
  },
}

local python_keys = {
  nodes = {
    expression_statement = true,
    parameters = "block",
    boolean_operator = "block",
    not_operator = "block",
    comparison_operator = "block",
    identifier = "block",
    call = "block",
    subscript = "block",
    with_clause = "block",
    pattern_list = "block",
  },
  options = {
    sibling_block_placement_behaviour = "before",
    placement_default_behaviour = "behind",
    language_format_behaviour = "full",
    print_function = "console.log",
  },
}

M.supported_languages = {
  lua = lua_keys,
  typescript = js_keys,
  typescriptreact = js_keys,
  javascript = js_keys,
  javascriptreact = js_keys,
  python = python_keys,
}

M.setup = function(opts)
  for k, v in pairs(opts) do
    --[[ can you spread tables in lua {...v}? ]]
    for inner_k, inner_v in pairs(v) do
      M.supported_languages[k]["options"][inner_k] = inner_v
    end
  end
end

local allowed_siblings = {
  block = true,
  statement_block = true,
}

local function getSibling(node)
  local type_text_node = node:type()
  local sibling = node:next_sibling()

  if sibling == nil then
    return node
  end

  while true do
    if allowed_siblings[sibling:type()] == nil then
      sibling = sibling:next_sibling()
      if sibling == nil then
        return node
      end
    else
      break
    end
  end

  local type_text_sibling = sibling:type()

  if M.language_keys["nodes"][type_text_node] == type_text_sibling then
    return sibling, M.language_keys["options"]["sibling_block_placement_behaviour"]
  else
    return node
  end
end

local function getParent(node, cnt)
  if cnt > 25 or node == nil then
    return nil
  end

  local type_text = node:type()

  local sibling_reference = M.language_keys["nodes"][type_text]
  if sibling_reference then
    if type(sibling_reference) == "boolean" then
      return node
    else
      local outnode, direction = getSibling(node)
      if direction == nil then
        return getParent(node:parent(), cnt + 1)
      else
        return outnode, direction
      end
    end
  end

  return getParent(node:parent(), cnt + 1)
end

M.get_gps = function()
  local status_gps_ok, gps = pcall(require, "nvim-navic")
  if not status_gps_ok then
    return ""
  end

  local status_ok, gps_location = pcall(gps.get_location, {})
  if not status_ok then
    return ""
  end

  if not gps.is_available() or gps_location == "error" then
    return ""
  end

  local breadcrumbs = "~"
  local sep = "%*"
  for i in string.gmatch(gps_location, "([^" .. sep .. "]+)") do
    local var = string.find(i, "NavicText")
    if var ~= nil then
      i = i:gsub("NavicText", "")
      i = i:gsub("%#", ""):gsub("%%", "")
      breadcrumbs = breadcrumbs .. i .. "->"
    end
  end
  return breadcrumbs:sub(1, -3)
end

M.get_text = function()
  local filename = vim.fn.expand("%:t")
  local current_text = vim.fn.expand("<cword>")
  local line_nr, _ = unpack(vim.api.nvim_win_get_cursor(0))

  local breadcrumbs = M.get_gps()

  return {
    filename = filename,
    line_nr = line_nr,
    breadcrumbs = breadcrumbs,
    current_text = current_text,
  }
end

M.get_text()

local function set_print_statement(node, print_text, placement, format) -- hello lll
  local start_row, start_idx, end_row, _ = node:range()

  local row_s
  local row_e
  local updated_text

  if placement == "inside" then
    row_s = start_row
    row_e = start_row + 1

    local current_text = vim.api.nvim_buf_get_lines(0, row_s, row_e, true)[1]
    updated_text = { current_text:sub(1, start_idx + 1) .. print_text .. ";" .. current_text:sub(start_idx + 2) }
  elseif placement == "before" then
    row_s = start_row
    row_e = start_row + 1

    local current_text = vim.api.nvim_buf_get_lines(0, row_s, row_e, true)[1]
    updated_text = current_text:sub(1, start_idx + 0) .. print_text .. ";"
    updated_text = { updated_text, current_text }
  elseif placement == "behind" then
    row_s = end_row
    row_e = end_row + 1

    local current_text = vim.api.nvim_buf_get_lines(0, row_s, row_e, true)[1]
    updated_text = { current_text .. ";" .. print_text }
  else
    row_s = end_row + 1
    row_e = end_row + 1
    updated_text = { print_text }
  end

  vim.api.nvim_buf_set_lines(0, row_s, row_e, true, updated_text)

  if format == "range" then
    vim.lsp.buf.format({ range = {
      ["start"] = { row_s + 0, 0 },
      ["end"] = { row_e + 2, -1 },
    } })
  elseif format == "full" then
    vim.lsp.buf.format({ async = true })
  end
end

M.insert_print = function()
  local ft = vim.bo.filetype
  M.language_keys = M.supported_languages[ft]
  if M.language_keys == nil then
    print("No language keys found for filetype: " .. ft)
    return
  end

  local node_at_cursor = ts_utils.get_node_at_cursor(0)
  local outp, placement = getParent(node_at_cursor, 0) --[[ or node_at_cursor ]]

  outp = outp or node_at_cursor
  placement = placement or M.language_keys["options"]["placement_default_behaviour"]
  local format = M.language_keys["options"]["language_format_behaviour"]
  local print_fn = M.language_keys["options"]["print_function"]
  local format_fn = M.language_keys["options"]["format_fn"]
  local t = M.get_text()
  local print_text = print_fn .. "(" .. format_fn(t) .. ")"

  set_print_statement(outp, print_text, placement, format)
end

return M

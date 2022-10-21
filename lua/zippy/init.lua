local M = {}

local ts_utils = require("nvim-treesitter.ts_utils")
local ts_locals = require("nvim-treesitter.locals")
local ts_query = require("nvim-treesitter.query")

-- in some helper module
function utils_Set(list)
	local set = {}
	for _, l in ipairs(list) do
		set[l] = true
	end
	return set
end

local lua_keys = {
	assignment_statement = true,
	function_call = true,
	parameters = "block",
	condition = "block",
}

--[[ local ts_statusline = require("nvim-treesitter.statusline") ]]
--[[  local f = require'nvim-treesitter'.statusline({}) ]]
--[[ print(f, "F") ]]

local usage_namespace = vim.api.nvim_create_namespace("nvim-treesitter-usages")

local cursor_pos = vim.api.nvim_win_get_cursor(0)

vim.api.nvim_set_hl(0, "ZIPPY", { default = true, bg = "#000000", fg = "#ffffff" })

local node_at_cursor = ts_utils.get_node_at_cursor(0)
local text = vim.treesitter.query.get_node_text(node_at_cursor, 0)
--[[ P(text) ]]

--[[ local parent = node_at_cursor:parent() ]]
--[[ local textParent = vim.treesitter.query.get_node_text(parent, 0) ]]
--[[ P(textParent) ]]

--[[ local parent = ts_utils.get_root_for_node(node_at_cursor) ]]
--[[ local textParent = vim.treesitter.query.get_node_text(parent, 0) ]]
--[[ P(textParent) ]]

--[[ ts_utils.highlight_node(node_at_cursor, 0, usage_namespace, "ZIPPY") ]]
--[[ print(node_at_cursor) ]]
--[[ print(ts_utils.get_node_range(node_at_cursor), "node_length") ]]
--[[ local range_of_node = ts_utils.get_node_range(node_at_cursor) ]]
--[[ P(range_of_node) ]]
--[[ ts_utils.goto_node(node_at_cursor, 1, 1) ]]
--[[ case variable ]]
--[[ start at itentifier and iterate until the next assignemt_statement, local decleration ]]

local function getSibling(node)
	local type_text_node = node:type()
	local sibling = node:next_sibling()

	if sibling == nil then
		return node
	end

	local type_text_sibling = sibling:type()

	if lua_keys[type_text_node] == type_text_sibling then
		return node:next_sibling()
	else
		return node
	end
end

local function getParent(node, cnt)
	if cnt > 5 or node == nil then
		return nil
	end

	local type_text = node:type()

	if lua_keys[type_text] then
		return node --[[ getSibling(node) ]]
	end

	return getParent(node:parent(), cnt + 1)
end

local function set_print_statement(node, text, placement, format)
	local start_row, _, end_row, _ = node:range()
	local row
	if placement == "above" then
		row = start_row
	else
		row = end_row + 1
	end
	vim.api.nvim_buf_set_lines(0, row, row, true, { text })
	if format then
		vim.lsp.buf.format({ range = {
			["start"] = { end_row + 2, 0 },
			["end"] = { end_row + 2, -1 },
		} })
	end
end

local test_variable, test2 = "test", "test2"

local outp = getParent(node_at_cursor, 0) or node_at_cursor
set_print_statement(outp, "print('HELLO')", "below", true)

--[[ local sibling1 = outp:next_sibling() ]]
--[[ print(sibling1:type(), "SIBLING1") ]]

--[[ local range_of_node = ts_utils.node_to_lsp_range(outp) ]]
--[[ P(range_of_node) ]]
--[[ ts_utils.goto_node(outp, 1, 1) ]]

--[[ 
-- break on:
--
-- assignment_statement
-- function_call
-- parameters
--
--]]

--[[ local res123 = ts_locals.iter_scope_tree(node_at_cursor, 0) ]]
--[[ for scope in res123 do ]]
--[[ 	print(scope, "scope") ]]
--[[ end ]]

--[[ print(node_at_cursor:parent():type(), "test_variable") ]]
--[[ print(node_at_cursor:parent():parent():type(), "test_variable") ]]

M.example = function(axxxxxxx, ababm, xccc)
	print(axxxxxxx)
	local blalalalalla = 10
	local function help(abc)
		print("HELLO")
	end
end

--[[ local r = ts_locals.get_scope_tree(node_at_cursor) ]]
--[[ local closest_scope = r[1] ]]
--[[ print(closest_scope, "closest_scope") ]]
--[[ print(node_at_cursor, "node_at_cursor") ]]

--[[ for k, v in ipairs(r) do ]]
--[[ 	print(k, v) ]]
--[[ 	P(v) ]]
--[[ 	local text_current = vim.treesitter.query.get_node_text(v, 0) ]]
--[[ 	print(text_current) ]]
--[[ 	ts_utils.highlight_node(v, 0, usage_namespace, "ZIPPY") ]]
--[[ 	break ]]
--[[ end ]]

--[[ ts_utils.highlight_node(current_scope, 0, usage_namespace, "ZIPPY") ]]
--[[ P(current_scope) ]]

--[[ local text_current = vim.treesitter.query.get_node_text(current_scope, 0) ]]
--[[ P(text_current) ]]

--[[ local res = ts_locals.get_scopes() ]]
--[[ P(res) ]]

return M

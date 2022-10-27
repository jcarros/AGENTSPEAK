# -*- coding: utf-8 -*-
#
# This file is part of the python-agentspeak interpreter.
# Copyright (C) 2016-2019 Niklas Fiekas <niklas.fiekas@tu-clausthal.de>.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function

import sys
import errno
import os.path

import agentspeak
import agentspeak.lexer
import agentspeak.util

from agentspeak import Trigger, GoalType, FormulaType, UnaryOp, BinaryOp


class AstBaseVisitor(object):
    def visit_literal(self, ast_literal):
        pass

    def visit_list(self, ast_list):
        pass

    def visit_linked_list(self, ast_linked_list):
        pass

    def visit_rule(self, ast_rule):
        pass

    def visit_goal(self, ast_goal):
        pass

    def visit_formula(self, ast_formula):
        pass

    def visit_const(self, ast_const):
        pass

    def visit_variable(self, ast_variable):
        pass

    def visit_unary_op(self, ast_unary_op):
        pass

    def visit_binary_op(self, ast_binary_op):
        pass

    def visit_plan(self, ast_plan):
        pass

    def visit_event(self, ast_event):
        pass

    def visit_body(self, ast_body):
        pass

    def visit_while(self, ast_while):
        pass

    def visit_for(self, ast_for):
        pass

    def visit_if_then_else(self, ast_if_then_else):
        pass

    def visit_agent(self, ast_agent):
        pass


class AstNode(object):
    def __init__(self):
        self.loc = None


class AstLiteral(AstNode):
    def __init__(self): 
        super(AstLiteral, self).__init__() 
        self.functor = None
        self.terms = []
        self.annotations = []

    def accept(self, visitor):
        return visitor.visit_literal(self)

    def signature(self):
        return "%s/%d" % (self.functor, len(self.terms))

    def __str__(self):
        builder = []
        builder.append(self.functor)
        if self.terms:
            builder.append("(")
            builder.append(", ".join(str(term) for term in self.terms))
            builder.append(")")
        if self.annotations:
            builder.append("[")
            builder.append(", ".join(str(term) for term in self.annotations))
            builder.append("]")
        return "".join(builder)


class AstList(AstNode):
    def __init__(self):
        super(AstList, self).__init__()
        self.terms = []

    def accept(self, visitor):
        return visitor.visit_list(self)

    def __str__(self):
        return "[%s]" % (", ".join(str(term) for term in self.terms), )


class AstLinkedList(AstNode):
    def __init__(self):
        super(AstLinkedList, self).__init__()
        self.head = None
        self.tail = None

    def accept(self, visitor):
        return visitor.visit_linked_list(self)

    def __str__(self):
        return "[%s|%s]" % (self.head, self.tail)


class AstRule(AstNode):
    def __init__(self):
        super(AstRule, self).__init__()
        self.head = None
        self.consequence = None

    def accept(self, visitor):
        return visitor.visit_rule(self)

    def __str__(self):
        return "%s :- %s" % (self.head, self.consequence)


class AstGoal(AstNode):
    def __init__(self):
        super(AstGoal, self).__init__()
        self.atom = None

    def accept(self, visitor):
        return visitor.visit_goal(self)

    def __str__(self):
        return "!%s" % (self.atom, )


class AstFormula(AstNode):
    def __init__(self):
        super(AstFormula, self).__init__()
        self.formula_type = None
        self.term = None

    def accept(self, visitor):
        return visitor.visit_formula(self)

    def __str__(self):
        return "%s%s" % (self.formula_type.value, str(self.term))


class AstConst(AstNode):
    def __init__(self):
        super(AstConst, self).__init__()
        self.value = None

    def accept(self, visitor):
        return visitor.visit_const(self)

    def __str__(self):
        return agentspeak.asl_repr(self.value)


class AstVariable(AstNode):
    def __init__(self):
        super(AstVariable, self).__init__()
        self.name = None
        self.proven_bound = False
        self.proven_unbound = False

    def accept(self, visitor):
        return visitor.visit_variable(self)

    def __str__(self):
        return self.name


class AstUnaryOp(AstNode):
    def __init__(self):
        super(AstUnaryOp, self).__init__()
        self.operator = None
        self.operand = None

    def accept(self, visitor):
        return visitor.visit_unary_op(self)

    def __str__(self):
        return "(%s %s)" % (self.operator.value.lexeme, self.operand)


class AstBinaryOp(AstNode):
    def __init__(self):
        super(AstBinaryOp, self).__init__()
        self.operator = None
        self.left = None
        self.right = None

    def accept(self, visitor):
        return visitor.visit_binary_op(self)

    def __str__(self):
        return "(%s %s %s)" % (self.left, self.operator.value.lexeme, self.right)


class AstPlan(AstNode):
    def __init__(self):
        super(AstPlan, self).__init__()
        self.annotations = []
        self.event = None
        self.context = None
        self.body = None

    def accept(self, visitor):
        return visitor.visit_plan(self)

    def signature(self):
        return self.event.signature()

    def __str__(self):
        builder = []

        for annotation in self.annotations:
            builder.append("@")
            builder.append(str(annotation))
            builder.append("\n")

        builder.append(str(self.event))

        if self.context:
            builder.append(" : ")
            builder.append(str(self.context))

        if self.body:
            builder.append(" <-\n")
            builder.append(agentspeak.util.indent(str(self.body)))

        return "".join(builder)


class AstEvent(AstNode):
    def __init__(self):
        super(AstEvent, self).__init__()
        self.trigger = None
        self.goal_type = None
        self.head = None

    def signature(self):
        return "%s%s%s" % (self.trigger.value, self.goal_type.value, self.head.signature())

    def accept(self, visitor):
        return visitor.visit_event(self)

    def __str__(self):
        builder = []
        builder.append(self.trigger.value)
        builder.append(self.goal_type.value)
        builder.append(str(self.head))
        return "".join(builder)


class AstBody(AstNode):
    def __init__(self):
        super(AstBody, self).__init__()
        self.formulas = []

    def accept(self, visitor):
        return visitor.visit_body(self)

    def __str__(self):
        builder = []

        first = True
        for formula in self.formulas:
            if not first:
                builder.append(";\n")
            first = False

            builder.append(str(formula))

        return "".join(builder)


class AstWhile(AstNode):
    def __init__(self):
        super(AstWhile, self).__init__()
        self.condition = None
        self.body = None

    def accept(self, visitor):
        return visitor.visit_while(self)

    def __str__(self):
        builder = []
        builder.append("while (")
        builder.append(str(self.condition))
        builder.append(") {\n")
        if self.body and self.body.formulas:
            builder.append(agentspeak.util.indent(str(self.body)))
            builder.append(";\n")
        builder.append("}")
        return "".join(builder)


class AstFor(AstNode):
    def __init__(self):
        super(AstFor, self).__init__()
        self.generator = None
        self.body = None

    def accept(self, visitor):
        return visitor.visit_for(self)

    def __str__(self):
        builder = []
        builder.append("for (")
        builder.append(str(self.generator))
        builder.append(") {\n")
        if self.body and self.body.formulas:
            builder.append(agentspeak.util.indent(str(self.body)))
            builder.append(";\n")
        builder.append("}")
        return "".join(builder)


class AstIfThenElse(AstNode):
    def __init__(self):
        super(AstIfThenElse, self).__init__()
        self.condition = None
        self.if_body = None
        self.else_body = None

    def accept(self, visitor):
        return visitor.visit_if_then_else(self)

    def __str__(self):
        builder = []
        builder.append("if (")
        builder.append(str(self.condition))
        builder.append(") {\n")
        if self.if_body and self.if_body.formulas:
            builder.append(agentspeak.util.indent(str(self.if_body)))
            builder.append(";\n")
        if self.else_body and self.else_body.formulas:
            builder.append("} else {\n")
            builder.append(agentspeak.util.indent(str(self.else_body)))
            builder.append(";\n")
        builder.append("}")
        return "".join(builder)


class AstAgent(AstNode):
    def __init__(self):
        super(AstAgent, self).__init__()
        self.rules = []
        self.beliefs = []
        self.goals = []
        self.plans = []

    def accept(self, visitor):
        return visitor.visit_agent(self)

    def __str__(self):
        builder = []
        for rule in self.rules:
            if builder:
                builder.append("\n")
            builder.append(str(rule))
            builder.append(".")

        if self.beliefs:
            if builder:
                builder.append("\n")
            for belief in self.beliefs:
                if builder:
                    builder.append("\n")
                builder.append(str(belief))
                builder.append(".")

        if self.goals:
            if builder:
                builder.append("\n")
            for goal in self.goals:
                if builder:
                    builder.append("\n")
                builder.append(str(goal))
                builder.append(".")

        for plan in self.plans:
            if builder:
                builder.append("\n\n")
            builder.append(str(plan))
            builder.append(".")

        return "".join(builder)


def parse_literal(tok, tokens, log): 
    if not tok.token.functor:  #If tok is not a functor
        raise log.error("expected functor, got '%s'", tok.lexeme, loc=tok.loc) #Raise an error

    literal = AstLiteral() # create a new literal
    literal.functor = tok.lexeme # set the functor
    literal.loc = tok.loc # set the location of the functor in the literal object

    tok = next(tokens) # get the next token

    if tok.lexeme == "(": # if the next token is a "("
        while True: # while true
            tok = next(tokens) # get the next token
            tok, term = parse_term(tok, tokens, log) # parse the term and get the next token
            literal.terms.append(term) # add the term to the literal

            if tok.lexeme == ")": # if the next token is a ")"
                tok = next(tokens) # get the next token
                break # break the loop
            elif tok.lexeme == ",": # if the next token is a ","
                continue # continue the loop
            else:
                raise log.error("expected ')' or another argument for the literal, got '%s'", # raise an error
                                tok.lexeme, loc=tok.loc, extra_locs=[literal.loc]) # with the location of the literal

    if tok.lexeme == "[": # if the next token is a "["
        while True: # while true
            tok = next(tokens) # get the next token
            tok, term = parse_term(tok, tokens, log) # parse the term and get the next token
            literal.annotations.append(term) # add the term to the annotations of the literal

            if tok.lexeme == "]": # if the next token is a "]"
                tok = next(tokens) # get the next token
                break
            elif tok.lexeme == ",": # if the next token is a ","
                continue
            else:
                raise log.error("expected ']' or another annotation, got '%s'", tok.lexeme,
                                loc=tok.loc, extra_locs=[literal.loc])

    return tok, literal # return the next token and the literal


def parse_list(tok, tokens, log):
    """
    This function parses a list of terms.

    :param tok: The current token.
    :param tokens: The token stream.
    :param log: The logger.
    
    if we have a list of terms, 
    we explore the list until the find the end of the list
    if we find a "|" and we have less or equal than 2 terms in the list
    we parse the term after the "|" and we return the list of terms and 
    the term after the "|" 
    """
    
    if tok.lexeme != "[": # if the next token is not a "["
        raise log.error("expected '[' for list, got '%s'", tok.lexeme, loc=tok.loc) # raise an error, this is not a list

    ast_list = AstList() # create a new AstList object
    ast_list.loc = tok.loc  # set the location of the list

    tok = next(tokens) # get the next token

    while tok.lexeme != "]": # while the next token is not a "]"
        tok, term = parse_and_expr(tok, tokens, log) # parse the term and get the next token
        ast_list.terms.append(term) # add the term to the list

        if tok.lexeme == "|": # if the next token is a "|"
            if len(ast_list.terms) > 1: # if the list has more than one term
                # [a, b, c | d] is not a valid list
                raise log.error("do not mix ',' and '|' in list notation",
                                loc=tok.loc, extra_locs=[ast_list.loc]) # raise an error because it is not a valid list
            # [a | b] is a valid list
            # This is a valid list with 5 terms: [a, b, c, d, e] 
            """If the list has only one term, then the next term is the tail of the list, 
            so we return the next token and the list"""
            ast_linked_list = AstLinkedList() # create a new AstLinkedList object
            ast_linked_list.loc = ast_list.loc # set the location of the linked list
            ast_linked_list.head = term # set the head of the linked list
            tok, ast_linked_list.tail = parse_linked_list_tail(tok, tokens, log) # parse the tail of the linked list and get the next token
            return tok, ast_linked_list # return the next token and the linked list
        elif tok.lexeme == ",": # if the next token is a ","
            tok = next(tokens) # get the next token
            continue
        elif tok.lexeme != "]": # if the next token is not a "]"
            raise log.error("expected ']' or another term for list, got '%s'", tok.lexeme,
                            loc=tok.loc, extra_locs=[ast_list.loc]) # raise an error because it is not a valid list notation

    return next(tokens), ast_list # return the next token and the list


def parse_linked_list_tail(tok, tokens, log):
    if tok.lexeme != "|": # if the next token is not a "|" raise an error
        # [a, b, c | d] is not a valid list
        raise log.error("expected '|' before linked list tail, got '%s'", tok.lexeme, loc=tok.loc)

    tok = next(tokens) # get the next token
    tok, term = parse_and_expr(tok, tokens, log) # parse the term and get the next token

    if tok.lexeme == "|": # if the next token is a "|" 
        ast_linked_list = AstLinkedList() # create a new AstLinkedList object
        ast_linked_list.loc = tok.loc # set the location of the linked list
        ast_linked_list.head = term # set the head of the linked list
        tok, ast_linked_list.tail = parse_linked_list_tail(tok, tokens, log) # parse the tail of the linked list and get the next token
        return tok, ast_linked_list # return the next token and the linked list
    elif tok.lexeme == "]":
        return next(tokens), term
    else:
        raise log.error("expected ']' or '|' followed by another term, got '%s'",
                        tok.lexeme, loc=tok.loc)


def parse_atom(tok, tokens, log):
    if tok.token.variable: # if the next token is a variable
        """
        If the next token is a variable, we create a new AstVariable 
        object and we set the location of the variable
        """
        variable = AstVariable() # create a new AstVariable object
        variable.name = tok.lexeme # set the name of the variable
        variable.loc = tok.loc # set the location of the variable
        return next(tokens), variable # return the next token and the variable
    elif tok.token.functor: # if the next token is a functor
        return parse_literal(tok, tokens, log) # parse the literal and return the next token and the literal
    elif tok.lexeme == "[": # if the next token is a "["
        return parse_list(tok, tokens, log) # parse the list and return the next token and the list
    elif tok.lexeme == "(": # if the next token is a "("
        tok = next(tokens) # get the next token
        tok, term = parse_term(tok, tokens, log) # parse the term and get the next token
        if tok.lexeme != ")": # if the next token is not a ")"
            raise log.error("expected ')' after term, got '%s'", tok.lexeme,
                            loc=tok.loc, extra_locs=[term.loc]) # raise an error because it is not a valid term
        return next(tokens), term # return the next token and the term
    elif tok.token.numeric: # if the next token is a numeric
        const = AstConst() # create a new AstConst object
        const.value = float(tok.lexeme) # set the value of the constant
        const.loc = tok.loc # set the location of the constant
        return next(tokens), const # return the next token and the constant
    elif tok.token.boolean is not None: # if the next token is a boolean
        const = AstConst() # create a new AstConst object
        const.value = tok.token.boolean # set the value of the constant
        const.loc = tok.loc # set the location of the constant
        return next(tokens), const # return the next token and the constant
    elif tok.token.string: # if the next token is a string
        const = AstConst() # create a new AstConst object
        const.value = agentspeak.parse_string(tok.lexeme) # set the value of the constant
        const.loc = tok.loc # set the location of the constant
        return next(tokens), const # return the next token and the constant
    else:
        raise log.error("expected term, got '%s'", tok.lexeme, loc=tok.loc) # raise an error because it is not a valid term


def parse_power(tok, tokens, log):
    tok, expr = parse_atom(tok, tokens, log) # parse the atom and get the next token
    while tok.lexeme == "**": # while the next token is a "**"
        op = AstBinaryOp() # create a new AstBinaryOp object
        op.left = expr # set the left operand of the binary operation
        op.operator = BinaryOp.op_pow # set the operator of the binary operation
        op.loc = tok.loc # set the location of the binary operation
        tok = next(tokens) # get the next token
        tok, op.right = parse_factor(tok, tokens, log) # parse the factor and get the next token
        expr = op  # set the expression to the binary operation

    return tok, expr # return the next token and the expression


def parse_factor(tok, tokens, log): # parse the factor
    if tok.token.unary_op: # if the next token is a unary operator
        op = AstUnaryOp() # create a new AstUnaryOp object
        op.operator = tok.token.unary_op # set the operator of the unary operation
        op.loc = tok.loc # set the location of the unary operation
        tok = next(tokens) # get the next token
        tok, op.operand = parse_factor(tok, tokens, log) # parse the factor and get the next token
        return tok, op # return the next token and the unary operation
    else: # if the next token is not a unary operator
        return parse_power(tok, tokens, log) # parse the power and return the next token and the expression


def parse_product(tok, tokens, log): # parse the product
    tok, expr = parse_factor(tok, tokens, log) # parse the factor and get the next token
    while tok.token.mult_op: # while the next token is a multiplication operator
        op = AstBinaryOp() # create a new AstBinaryOp object
        op.left = expr # set the left operand of the binary operation
        op.operator = tok.token.mult_op # set the operator of the binary operation
        op.loc = tok.loc # set the location of the binary operation
        tok = next(tokens) # get the next token
        tok, op.right = parse_factor(tok, tokens, log) # parse the factor and get the next token
        expr = op # set the expression to the binary operation

    return tok, expr # return the next token and the expression


def parse_arith_expr(tok, tokens, log):
    tok, expr = parse_product(tok, tokens, log) # parse the product and get the next token
    while tok.token.add_op: # while the next token is an addition operator
        op = AstBinaryOp() # create a new AstBinaryOp object
        op.left = expr # set the left operand of the binary operation
        op.operator = tok.token.add_op # set the operator of the binary operation
        op.loc = tok.loc # set the location of the binary operation
        tok = next(tokens) # get the next token
        tok, op.right = parse_product(tok, tokens, log) # parse the product and get the next token
        expr = op # set the expression to the binary operation

    return tok, expr # return the next token and the expression


def parse_comparison(tok, tokens, log):
    comparisons = None
    tok, rightmost = parse_arith_expr(tok, tokens, log)
    while tok.token.comp_op:
        op = AstBinaryOp()
        op.left = rightmost
        op.operator = tok.token.comp_op
        op.loc = tok.loc
        tok = next(tokens)
        tok, rightmost = parse_arith_expr(tok, tokens, log)
        op.right = rightmost

        if comparisons is None:
            comparisons = op
        else:
            op_and = AstBinaryOp()
            op_and.left = comparisons
            op_and.operator = BinaryOp.op_and
            op_and.loc = op.loc
            op_and.right = op
            comparisons = op_and

    return tok, comparisons or rightmost


def parse_not_expr(tok, tokens, log):
    if tok.lexeme == "not":
        op = AstUnaryOp()
        op.operator = UnaryOp.op_not
        op.loc = tok.loc
        tok = next(tokens)
        tok, op.operand = parse_not_expr(tok, tokens, log)
        return tok, op

    return parse_comparison(tok, tokens, log)


def parse_and_expr(tok, tokens, log):
    tok, expr = parse_not_expr(tok, tokens, log)
    while tok.lexeme == "&":
        bin_op = AstBinaryOp()
        bin_op.left = expr
        bin_op.operator = BinaryOp.op_and
        bin_op.loc = tok.loc

        tok = next(tokens)
        tok, bin_op.right = parse_not_expr(tok, tokens, log)
        expr = bin_op

    return tok, expr


def parse_term(tok, tokens, log):
    tok, expr = parse_and_expr(tok, tokens, log)
    while tok.lexeme == "|":
        bin_op = AstBinaryOp()
        bin_op.left = expr
        bin_op.operator = BinaryOp.op_or
        bin_op.loc = tok.loc

        tok = next(tokens)
        tok, bin_op.right = parse_and_expr(tok, tokens, log)
        expr = bin_op

    return tok, expr


def parse_rule_or_belief(tok, tokens, log):
    if "." in tok.lexeme: 
        log.warning("found '.' in assertion. should this have been an action?", loc=tok.loc)

    tok, belief_atom = parse_literal(tok, tokens, log) # parse the literal and get the next token

    if tok.lexeme == ":-": # if the next token is a ":-" (i.e. if the next token is a rule)
        # A rule with head and body.
        rule = AstRule() # create a new AstRule object
        rule.head = belief_atom # set the head of the rule to the belief atom
        rule.loc = tok.loc # set the location of the rule

        tok = next(tokens) # get the next token
        tok, rule.consequence = parse_term(tok, tokens, log) # parse the term and get the next token
        return tok, rule # return the next token and the rule
    else:
        # Just the belief atom.
        return tok, belief_atom # return the next token and the belief atom


def parse_initial_goal(tok, tokens, log):
    if tok.lexeme != "!": # if the next token is not a "!"
        raise log.error("expected '!' for initial goal, got '%s'", tok.lexeme, loc=tok.loc) # raise an error

    goal = AstGoal() # create a new AstGoal object
    goal.loc = tok.loc # set the location of the goal

    tok = next(tokens) # get the next token
    tok, goal.atom = parse_literal(tok, tokens, log) # parse the literal and get the next token
    return tok, goal # return the next token and the goal


def parse_body(tok, tokens, log):
    body = AstBody() # create a new AstBody object
    body.loc = tok.loc # set the location of the body

    while tok.lexeme != "}": # while the next token is not a "}"
        tok, formula = parse_body_formula(tok, tokens, log) # parse the body formula and get the next token
        body.formulas.append(formula) # add the body formula to the body

        if tok.lexeme == ";": # if the next token is a ";"
            tok = next(tokens) # get the next token
            if tok.lexeme == "}": # if the next token is a "}"
                break
            continue
        elif tok.lexeme == "}":
            break
        elif isinstance(formula, AstFormula): # if the body formula is an AstFormula
            raise log.error("expected ';' or '}' after formula, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[formula.loc])
        else:
            # Block in brackets has ended.
            pass

    return tok, body


def parse_while(tok, tokens, log):
    if tok.lexeme != "while": # if the next token is not a "while"
        raise log.error("expected 'while', got '%s'", tok.lexeme, loc=tok.loc) # raise an error

    while_node = AstWhile() # create a new AstWhile object
    while_node.loc = tok.loc   # set the location of the while node
    tok = next(tokens) # get the next token

    if tok.lexeme != "(": # if the next token is not a "("
        raise log.error("expected '(' after while, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[while_node.loc])
    tok = next(tokens) # get the next token

    tok, while_node.condition = parse_term(tok, tokens, log) # parse the term and get the next token

    if tok.lexeme != ")": # if the next token is not a ")"
        raise log.error("expected ')' after while condition, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[while_node.loc, while_node.condition.loc])
    tok = next(tokens) # get the next token

    if tok.lexeme != "{": # if the next token is not a "{"
        raise log.error("expected '{' after while head, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[while_node.loc])
    tok = next(tokens) # get the next token

    tok, while_node.body = parse_body(tok, tokens, log) # parse the body and get the next token

    if tok.lexeme != "}": # if the next token is not a "}"
        raise log.error("expected '}' after while body, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[while_node.loc])
    tok = next(tokens) # get the next token

    return tok, while_node # return the next token and the while node


def parse_for(tok, tokens, log):
    for_node = AstFor() # create a new AstFor object
    for_node.loc = tok.loc # set the location of the for node

    if tok.lexeme != "for": # if the next token is not a "for"
        raise log.error("expected 'for', got '%s'", tok.lexeme, loc=tok.loc) # raise an error
    tok = next(tokens) # get the next token

    if tok.lexeme != "(": # if the next token is not a "("
        raise log.error("expected '(' after for, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[for_node.loc])
    tok = next(tokens) # get the next token

    tok, for_node.generator = parse_term(tok, tokens, log) # parse the term and get the next token

    if tok.lexeme != ")": # if the next token is not a ")"
        raise log.error("expected ')' after for generator, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[for_node.loc, for_node.generator.loc])
    tok = next(tokens) # get the next token

    if tok.lexeme != "{": # if the next token is not a "{"
        raise log.error("expected '{' after for head, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[for_node.loc])
    tok = next(tokens) # get the next token

    tok, for_node.body = parse_body(tok, tokens, log) # parse the body and get the next token

    if tok.lexeme != "}": # if the next token is not a "}"
        raise log.error("expected '}' after for body, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[for_node.loc])
    tok = next(tokens) # get the next token

    return tok, for_node # return the next token and the for node


def parse_if_then_else(tok, tokens, log):
    if_then_else = AstIfThenElse() # create a new AstIfThenElse object
    if_then_else.loc = tok.loc # set the location of the if then else

    if tok.lexeme != "if": # if the next token is not a "if"
        raise log.error("expected 'if', got '%s'", tok.lexeme, loc=tok.loc)
    tok = next(tokens) # get the next token

    if tok.lexeme != "(": # if the next token is not a "("
        raise log.error("expected '(' after if, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[if_then_else.loc])
    tok = next(tokens) # get the next token

    tok, if_then_else.condition = parse_term(tok, tokens, log) # parse the term and get the next token

    if tok.lexeme != ")": # if the next token is not a ")"
        raise log.error("expected ')' after if condition, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[if_then_else.loc, if_then_else.condition.loc])
    tok = next(tokens) # get the next token

    if tok.lexeme != "{": # if the next token is not a "{"
        raise log.error("expected '{' after if head, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[if_then_else.loc])
    tok = next(tokens) # get the next token

    tok, if_then_else.if_body = parse_body(tok, tokens, log) # parse the body and get the next token

    if tok.lexeme != "}": # if the next token is not a "}"
        raise log.error("expected '}' after if body, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[if_then_else.loc])
    tok = next(tokens) # get the next token

    if tok.lexeme == "else": # if the next token is a "else"
        tok_else = tok # set the else token
        tok = next(tokens) # get the next token

        if tok.lexeme != "{": # if the next token is not a "{"
            raise log.error("expected '{' after else, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[if_then_else.loc, tok_else.loc])
        tok = next(tokens) # get the next token

        tok, if_then_else.else_body = parse_body(tok, tokens, log) # parse the body and get the next token

        if tok.lexeme != "}": # if the next token is not a "}"
            raise log.error("expected '}' after else body, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[if_then_else.loc, tok_else.loc])
        tok = next(tokens) # get the next token

    return tok, if_then_else # return the next token and the if then else node


def parse_body_formula(tok, tokens, log):
    if tok.lexeme == "while": # if the next token is a "while"
        return parse_while(tok, tokens, log)
    elif tok.lexeme == "if": # if the next token is a "if"
        return parse_if_then_else(tok, tokens, log)
    elif tok.lexeme == "for": # if the next token is a "for"
        return parse_for(tok, tokens, log)
    else:
        formula = AstFormula() # create a new AstFormula object
        formula.loc = tok.loc # set the location of the formula

        if tok.token.formula_type: # if the next token is a formula type
            formula.formula_type = tok.token.formula_type  # set the formula type
            tok = next(tokens) # get the next token
        else:
            formula.formula_type = FormulaType.term # set the formula type to term

        tok, formula.term = parse_term(tok, tokens, log) # parse the term and get the next token
        return tok, formula # return the next token and the formula node


def parse_plan_body(tok, tokens, log):
    body = AstBody() # create a new AstBody object

    while True:
        tok, formula = parse_body_formula(tok, tokens, log) # parse the body formula and get the next token
        body.formulas.append(formula) # append the formula to the body

        if tok.lexeme == ";": # if the next token is a ";"
            tok = next(tokens) # get the next token
            continue
        elif tok.lexeme == ".": # if the next token is a "."
            break
        elif isinstance(formula, AstFormula): # if the formula is a AstFormula
            raise log.error("expected ';' or '.' after formula, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[formula.loc])
        else:
            # Block in brackets has ended.
            pass

    return tok, body # return the next token and the body node


def parse_event(tok, tokens, log):
    event = AstEvent() # create a new AstEvent object

    if not tok.token.trigger: # if the next token is not a trigger
        raise log.error("expected plan trigger, got '%s'", tok.lexeme, loc=tok.loc)
    event.loc = tok.loc # set the location of the event
    event.trigger = tok.token.trigger # set the trigger of the event
    tok = next(tokens) # get the next token

    if tok.token.goal_type: # if the next token is a goal type
        event.goal_type = tok.token.goal_type # set the goal type of the event
        tok = next(tokens) # get the next token
    else:
        event.goal_type = GoalType.belief # set the goal type to belief

    tok, event.head = parse_literal(tok, tokens, log) # parse the literal and get the next token
    return tok, event # return the next token and the event node


def parse_plan(tok, tokens, log):
    plan = AstPlan() # create a new AstPlan object

    while tok.lexeme == "@": # while the next token is a "@"
        tok = next(tokens) # get the next token
        tok, annotation = parse_literal(tok, tokens, log) # parse the literal and get the next token
        plan.annotations.append(annotation) # append the annotation to the plan

    tok, event = parse_event(tok, tokens, log) # parse the event and get the next token
    plan.event = event # set the event of the plan
    plan.loc = event.loc # set the location of the plan

    if tok.lexeme == ":": # if the next token is a ":"
        tok = next(tokens) # get the next token
        tok, plan.context = parse_term(tok, tokens, log) # parse the term and get the next token

    if tok.lexeme == "<-": # if the next token is a "<-"
        body_loc = tok.loc # set the location of the body
        tok = next(tokens) # get the next token

        tok, plan.body = parse_plan_body(tok, tokens, log) # parse the body and get the next token
        plan.body.loc = body_loc # set the location of the body

    return tok, plan # return the next token and the plan node


def parse_agent(filename, tokens, log, included_files, directive=None):
    included_files = included_files | frozenset([os.path.normpath(filename)]) # add the filename to the included files
    agent = AstAgent() # create a new AstAgent object
    last_plan = None # set the last plan to None

    while True:
        try:
            tok = next(tokens) # get the next token
        except StopIteration: 
            if directive:
                # TODO: Where was the directive started?
                raise log.error("end of file, but did not close directive '%s'", directive)
            return validate(agent, log) # return the agent
        """
        Cambios Jordi 25/10/20022
        """
        if tok.lexeme == "{": # if the next token is a "{"
            tok = next(tokens) # get the next token
            if tok.lexeme == "include": # if the next token is a "include"
                include_loc = tok.loc # set the location of the include
                tok = next(tokens) # get the next token
                if tok.lexeme != "(": # if the next token is not a "(" 
                    raise log.error("expected '(' after include, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[include_loc])
                tok = next(tokens) # get the next token
                if not tok.token.string: # if the next token is not a string
                    raise log.error("expected filename to include, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[include_loc])
                include = agentspeak.parse_string(tok.lexeme) # parse the string
                tok = next(tokens) # get the next token
                if tok.lexeme != ")": # if the next token is not a ")"
                    raise log.error("expected ')' after include filename, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[include_loc])
                tok = next(tokens) # get the next token
                if tok.lexeme != "}": # if the next token is not a "}"
                    raise log.error("expected '}' to close include directive, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[include_loc])

                # Resolve included path.
                include = os.path.join(os.path.dirname(filename), include) # get the path of the include

                # Parse included file.
                if include in included_files: # if the include is already included
                    log.error("infinite recursive include: '%s'", include, loc=include_loc)
                else:
                    try:
                        included_file = open(include) # open the include    
                    except IOError as err: # if the include does not exist
                        if err.errno == errno.ENOENT: # if the error is that the file does not exist
                            log.error("include file not found: '%s'", include, loc=include_loc) # log the error
                        else:
                            raise
                    else:
                        included_tokens = agentspeak.lexer.TokenStream(included_file, 1) # create a new token stream
                        included_agent = parse(include, included_tokens, log, included_files) # parse the include
                        agent.beliefs += included_agent.beliefs # add the beliefs of the included agent to the agent
                        agent.rules += included_agent.rules # add the rules of the included agent to the agent
                        agent.goals += included_agent.goals # add the goals of the included agent to the agent
                        agent.plans += included_agent.plans # add the plans of the included agent to the agent
                        included_file.close() # close the include
            elif tok.lexeme == "begin": # if the next token is a "begin"
                begin_loc = tok.loc # set the location of the begin
                tok = next(tokens) # get the next token
                tok, sub_directive = parse_literal(tok, tokens, log) # parse the literal and get the next token
                if tok.lexeme != "}": # if the next token is not a "}"
                    raise log.error("expected '}' after begin, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[begin_loc])
                log.warning("directives are ignored as of yet", loc=sub_directive.loc)
                sub_agent = parse(filename, tokens, log, included_files, sub_directive) # parse the sub agent
                agent.beliefs += sub_agent.beliefs # add the beliefs of the sub agent to the agent
                agent.rules += sub_agent.rules # add the rules of the sub agent to the agent
                agent.goals += sub_agent.goals # add the goals of the sub agent to the agent
                agent.plans += sub_agent.plans # add the plans of the sub agent to the agent
            elif tok.lexeme == "end": # if the next token is a "end"
                end_loc = tok.loc # set the location of the end
                tok = next(tokens) # get the next token
                if tok.lexeme != "}": # if the next token is not a "}"
                    raise log.error("expected '}' after end, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[end_loc])
                if not directive: # if there is no directive
                    log.error("unexpected end", loc=end_loc)
                else:
                    return validate(agent, log) # return the agent
            else:
                raise log.error("expected 'include', or 'begin' or 'end' after '{', got '%s'", tok.lexeme, loc=tok.loc)
        elif tok.token.functor: # if the next token is a functor
            if last_plan is not None: # if there is a last plan
                log.warning("assertion after plan. should this have been part of '%s'?", last_plan.signature(), loc=tok.loc)
            tok, ast_node = parse_rule_or_belief(tok, tokens, log) # parse the rule or belief and get the next token
            if isinstance(ast_node, AstRule): # if the ast node is a rule
                if tok.lexeme != ".": # if the next token is not a "."
                    log.info("missing '.' after this rule", loc=ast_node.loc)
                    raise log.error("expected '.' after rule, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[ast_node.loc])
                agent.rules.append(ast_node) # add the rule to the agent
            else:
                if tok.lexeme != ".": # if the next token is not a "."
                    print("Here is the error")
                    log.info("missing '.' after this belief", loc=ast_node.loc)
                    raise log.error("expected '.' after belief, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[ast_node.loc])
                agent.beliefs.append(ast_node) # add the belief to the agent
        elif tok.lexeme == "hola":
            print("hola")
        elif tok.lexeme == "!": # if the next token is a "!"
            tok, ast_node = parse_initial_goal(tok, tokens, log) # parse the initial goal and get the next token
            if tok.lexeme != ".": # if the next token is not a "."
                log.info("missing '.' after this goal", loc=ast_node.loc)
                raise log.error("expected '.' after initial goal, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[ast_node.loc])
            agent.goals.append(ast_node) # add the goal to the agent
        elif tok.lexeme in ["@", "+", "-"]: # if the next token is a "@", "+" or "-"
            tok, last_plan = parse_plan(tok, tokens, log) # parse the plan and get the next token
            if tok.lexeme != ".": # if the next token is not a "."
                log.info("missing '.' after this plan", loc=last_plan.loc)
                raise log.error("expected '.' after plan, got '%s'", tok.lexeme, loc=tok.loc, extra_locs=[last_plan.loc])
            agent.plans.append(last_plan) # add the plan to the agent
        else:
            log.error("unexpected token: '%s'", tok.lexeme, loc=tok.loc) # log the error


class FindVariablesVisitor(object):
    def visit_literal(self, ast_literal):
        for term in ast_literal.terms: # for each term in the literal
            for var in term.accept(self): # for each variable in the term
                yield var # yield the variable

        for annotation in ast_literal.annotations: # for each annotation in the literal
            for var in annotation.accept(self): # for each variable in the annotation
                yield var # yield the variable

    def visit_list(self, ast_list): 
        for term in ast_list.terms: # for each term in the list
            for var in term.accept(self): # for each variable in the term
                yield var # yield the variable

    def visit_const(self, ast_const):
        return
        yield

    def visit_variable(self, ast_variable):
        yield ast_variable

    def visit_unary_op(self, unary_op):
        for var in unary_op.operand.accept(self): # for each variable in the operand
            yield

    def visit_binary_op(self, binary_op):
        for var in binary_op.left.accept(self): # for each variable in the left
            yield

        for var in binary_op.right.accept(self): # for each variable in the right
            yield


class FindOpVisitor(object):
    def visit_literal(self, ast_literal):
        for term in ast_literal.terms: # for each term in the literal
            for op in term.accept(self): # for each operator in the term
                yield op # yield the operator

        for annotation in ast_literal.annotations: # for each annotation in the literal
            for op in annotation.accept(self): # for each operator in the annotation
                yield op # yield the operator

    def visit_list(self, ast_list):
        for term in ast_list.terms: # for each term in the list
            for op in term.accept(self): # for each operator in the term
                yield op # yield the operator

    def visit_linked_list(self, ast_linked_list):
        for op in ast_linked_list.head.accept(self): # for each operator in the head
            yield op # yield the operator
        for op in ast_linked_list.tail.accept(self): # for each operator in the tail
            yield op # yield the operator

    def visit_const(self, ast_const):
        return
        yield

    def visit_variable(self, ast_variable):
        return
        yield

    def visit_unary_op(self, unary_op):
        yield unary_op # yield the operator

    def visit_binary_op(self, binary_op):
        yield binary_op # yield the operator


class NumericFoldVisitor(object):
    def __init__(self, log):
        self.log = log # set the log

    def visit_binary_op(self, ast_binary_op):
        if ast_binary_op.operator.value.numeric_op: # if the operator is a numeric operator
            left = ast_binary_op.left.accept(self) # get the left
            right = ast_binary_op.right.accept(self) # get the right
            if (isinstance(left, AstConst) and agentspeak.is_number(left.value) and
                    isinstance(right, AstConst) and agentspeak.is_number(right.value)): # if the left and right are both constants and numbers
                try:
                    const = AstConst() # create a constant
                    const.loc = ast_binary_op.loc # set the location
                    const.value = ast_binary_op.operator.value.func(left.value, right.value) # set the value
                    return const # return the constant
                except ZeroDivisionError as err: # if there is a zero division error
                    self.log.error("%s", err, loc=ast_binary_op.loc, extra_locs=[left.loc, right.loc])
            else: # if the left or right is not a constant or a number
                ast_binary_op.left = left # set the left
                ast_binary_op.right = right # set the right
        else: # if the operator is not a numeric operator
            self.log.error("unexpected operator '%s' in numeric context",
                           ast_binary_op.operator.value.lexeme,
                           loc=ast_binary_op.loc,
                           extra_locs=[ast_binary_op.left.loc, ast_binary_op.right.loc])

        return ast_binary_op # return the binary operator

    def visit_unary_op(self, ast_unary_op):
        if ast_unary_op.operator.value.numeric_op: # if the operator is a numeric operator
            folded = ast_unary_op.operand.accept(self) # get the folded operand
            if isinstance(folded, AstConst) and agentspeak.is_number(folded.value): # if the folded operand is a constant and a number
                const = AstConst() # create a constant
                const.loc = ast_unary_op.loc # set the location
                const.value = ast_unary_op.operator.value.func(folded.value) # set the value
                return const # return the constant
            else: # if the folded operand is not a constant or a number
                ast_unary_op.operand = folded # set the operand
        else: # if the operator is not a numeric operator
            self.log.error("unexpected operator '%s' in numeric context",
                           ast_unary_op.operator.value.lexeme,
                           loc=ast_unary_op.loc,
                           extra_locs=[ast_unary_op.operand.loc])

        return ast_unary_op # return the unary operator

    def visit_variable(self, ast_variable):
        return ast_variable # return the variable

    def visit_const(self, ast_const):
        if ast_const.value is True or ast_const.value is False: # if the constant is a boolean
            self.log.error("boolean in numeric context", loc=ast_const.loc)
        elif isinstance(ast_const.value, str): # if the constant is a string
            self.log.error("string in numeric context", loc=ast_const.loc)

        return ast_const # return the constant

    def visit_literal(self, ast_literal):
        self.log.error("did not expect literal in numeric context", loc=ast_literal.loc)
        return ast_literal # return the literal

    def visit_list(self, ast_list):
        self.log.error("did not expect list in numeric context", loc=ast_list.loc)
        return ast_list # return the list


class BooleanFoldVisitor(object):
    def __init__(self, log):
        self.log = log # set the log

    def visit_binary_op(self, ast_binary_op):
        if ast_binary_op.operator.value.boolean_op: # if the operator is a boolean operator
            left = ast_binary_op.left.accept(self) # get the left
            right = ast_binary_op.right.accept(self) # get the right
            if (isinstance(left, AstConst) and isinstance(left.value, bool) and
                    isinstance(right, AstConst) and isinstance(right.value, bool)): # if the left and right are both constants and booleans
                const = AstConst()  # create a constant
                const.loc = ast_binary_op.loc # set the location
                const.value = ast_binary_op.operator.value.func(left.value, right.value) # set the value
                return const # return the constant
            else: # if the left or right is not a constant or a boolean
                ast_binary_op.left = left # set the left
                ast_binary_op.right = right # set the right
        elif ast_binary_op.operator.value.comp_op: # if the operator is a comparison operator
            left = ast_binary_op.left.accept(TermFoldVisitor(self.log)) # get the left 
            right = ast_binary_op.right.accept(TermFoldVisitor(self.log)) # get the right
            if isinstance(left, AstConst) and isinstance(right, AstConst): # if the left and right are both constants
                const = AstConst() # create a constant
                const.loc = ast_binary_op.loc # set the location
                if agentspeak.is_number(left.value) and agentspeak.is_number(right.value): # if the left and right are both numbers
                    const.value = ast_binary_op.operator.value.func(left.value, right.value) # set the value
                    return const
                elif isinstance(left.value, bool) and isinstance(right.value, bool): # if the left and right are both booleans
                    const.value = ast_binary_op.operator.value.func(left.value, right.value) # set the value
                    return const
                elif isinstance(left.value, str) and isinstance(right.value, str): # if the left and right are both strings
                    const.value = ast_binary_op.operator.value.func(left.value, right.value) # set the value
                    return const

            ast_binary_op.left = left # set the left
            ast_binary_op.right = right # set the right
            return ast_binary_op # return the binary operator
        else: # if the operator is not a boolean or comparison operator
            self.log.error("unexpected operator '%s' in boolean context",
                           ast_binary_op.operator,
                           loc=ast_binary_op.loc,
                           extra_locs=[ast_binary_op.left.loc, ast_binary_op.right.loc])

        return ast_binary_op # return the binary operator

    def visit_unary_op(self, ast_unary_op):
        if ast_unary_op.operator.value.boolean_op: # if the operator is a boolean operator
            folded = ast_unary_op.operand.accept(self) # get the folded operand
            if isinstance(folded, AstConst) and isinstance(folded.value, bool): # if the folded operand is a constant and a boolean
                const = AstConst() # create a constant
                const.loc = ast_unary_op.loc # set the location
                const.value = ast_unary_op.operator.value.func(folded.value) # set the value
                return const # return the constant
            else: # if the folded operand is not a constant or a boolean
                ast_unary_op.operand = folded # set the operand
        else: # if the operator is not a boolean operator
            self.log.error("unexpected operator '%s' in boolean context",
                           ast_unary_op.operator,
                           loc=ast_unary_op.loc,
                           extra_locs=[ast_unary_op.operand.loc])

        return ast_unary_op # return the unary operator

    def visit_variable(self, ast_variable):
        return ast_variable # return the variable

    def visit_const(self, ast_const):
        if isinstance(ast_const.value, str): # if the constant is a string
            self.log.error("string in boolean context", loc=ast_const.loc)
        elif agentspeak.is_number(ast_const.value): # if the constant is a number
            self.log.error("number '%s' in boolean context", ast_const.value, loc=ast_const.loc)

        return ast_const # return the constant

    def visit_literal(self, ast_literal):
        self.log.error("literal in boolean context", loc=ast_literal.loc)
        return ast_literal # return the literal

    def visit_list(self, ast_list):
        self.log.error("did not expect list in boolean context", loc=ast_list.loc)


class TermFoldVisitor(object):
    def __init__(self, log):
        self.log = log # set the log

    def visit_binary_op(self, ast_binary_op):
        if ast_binary_op.operator.value.numeric_op: # if the operator is a numeric operator
            return ast_binary_op.accept(NumericFoldVisitor(self.log)) # return the numeric fold
        else: # if the operator is not a numeric operator
            return ast_binary_op.accept(BooleanFoldVisitor(self.log))  # return the boolean fold

    def visit_unary_op(self, ast_unary_op):
        if ast_unary_op.operator.value.boolean_op: # if the operator is a boolean operator
            return ast_unary_op.accept(BooleanFoldVisitor(self.log)) # return the boolean fold
        else: # if the operator is not a boolean operator
            return ast_unary_op.accept(NumericFoldVisitor(self.log)) # return the numeric fold

    def visit_variable(self, ast_variable):
        return ast_variable # return the variable

    def visit_const(self, ast_const):
        return ast_const # return the constant

    def visit_literal(self, ast_literal):
        ast_literal.terms = [term.accept(self) for term in ast_literal.terms] # get the terms
        ast_literal.annotations = [annotation.accept(self) for annotation in ast_literal.annotations] # get the annotations
        return ast_literal # return the literal

    def visit_list(self, ast_list):
        ast_list.terms = [term.accept(self) for term in ast_list.terms] # get the terms
        return ast_list # return the list

    def visit_linked_list(self, ast_linked_list):
        ast_linked_list.head = ast_linked_list.head.accept(self) # get the head
        ast_linked_list.tail = ast_linked_list.tail.accept(self) # get the tail
        return ast_linked_list # return the linked list


class LogicalFoldVisitor(BooleanFoldVisitor):

    def visit_binary_op(self, ast_binary_op):
        if ast_binary_op.operator.value.query_op and not ast_binary_op.operator.value.boolean_op: # if the operator is a query operator and not a boolean operator
            ast_binary_op.left = ast_binary_op.left.accept(TermFoldVisitor(self.log)) # get the left
            ast_binary_op.right = ast_binary_op.right.accept(TermFoldVisitor(self.log)) # get the right
            return ast_binary_op # return the binary operator
        else: # if the operator is not a query operator or a boolean operator
            return super(LogicalFoldVisitor, self).visit_binary_op(ast_binary_op) # return the binary operator

    def visit_literal(self, ast_literal):
        return ast_literal.accept(TermFoldVisitor(self.log)) # return the literal


class ConstFoldVisitor(object):
    def __init__(self, log):
        self.log = log # set the log

    def visit_binary_op(self, ast_binary_op):
        return ast_binary_op.accept(TermFoldVisitor(self.log)) # return the binary operator

    def visit_unary_op(self, ast_unary_op):
        return ast_unary_op.accept(TermFoldVisitor(self.log)) # return the unary operator

    def visit_agent(self, ast_agent):
        ast_agent.rules = [rule.accept(self) for rule in ast_agent.rules] # get the rules
        ast_agent.beliefs = [belief.accept(self) for belief in ast_agent.beliefs] # get the beliefs
        ast_agent.goals = [goal.accept(self) for goal in ast_agent.goals] # get the goals
        ast_agent.plans = [plan.accept(self) for plan in ast_agent.plans] # get the plans
        return ast_agent # return the agent

    def visit_if_then_else(self, ast_if_then_else):
        ast_if_then_else.condition = ast_if_then_else.condition.accept(LogicalFoldVisitor(self.log)) # get the condition
        ast_if_then_else.if_body = ast_if_then_else.if_body.accept(self) # get the if body
        ast_if_then_else.else_body = ast_if_then_else.else_body.accept(self) if ast_if_then_else.else_body else None # get the else body
        return ast_if_then_else # return the if then else

    def visit_for(self, ast_for):
        ast_for.generator = ast_for.generator.accept(LogicalFoldVisitor(self.log)) # get the generator
        ast_for.body = ast_for.body.accept(self) # get the body
        return ast_for # return the for

    def visit_while(self, ast_while):
        ast_while.condition = ast_while.condition.accept(LogicalFoldVisitor(self.log)) # get the condition
        ast_while.body = ast_while.body.accept(self) # get the body
        return ast_while # return the while

    def visit_body(self, ast_body):
        ast_body.formulas = [formula.accept(self) for formula in ast_body.formulas] # get the formulas
        return ast_body # return the body
  
    def visit_event(self, ast_event):
        ast_event.head = ast_event.head.accept(TermFoldVisitor(self.log)) # get the head
        return ast_event # return the event

    def visit_plan(self, ast_plan):
        ast_plan.annotations = [annotation.accept(TermFoldVisitor(self.log)) for annotation in ast_plan.annotations] # get the annotations
        ast_plan.event = ast_plan.event.accept(self) # get the event
        ast_plan.context = ast_plan.context.accept(LogicalFoldVisitor(self.log)) if ast_plan.context else None # get the context
        ast_plan.body = ast_plan.body.accept(self) if ast_plan.body else None # get the body
        return ast_plan # return the plan

    def visit_variable(self, ast_variable):
        return ast_variable # return the variable

    def visit_const(self, ast_const):
        return ast_const # return the constant

    def visit_formula(self, ast_formula):
        if ast_formula.formula_type == FormulaType.term: # if the formula is a term
            ast_formula.term = ast_formula.term.accept(LogicalFoldVisitor(self.log)) # get the term
        else: # if the formula is not a term
            if isinstance(ast_formula.term, (AstLiteral, AstVariable)): # if the term is a literal or a variable
                ast_formula.term = ast_formula.term.accept(TermFoldVisitor(self.log)) # get the term
            else: # if the term is not a literal or a variable
                self.log.error("expected literal or variable after '%s'", ast_formula.formula_type, loc=ast_formula.loc, extra_locs=[ast_formula.term.loc]) 

        return ast_formula # return the formula

    def visit_goal(self, ast_goal):
        ast_goal.atom = ast_goal.atom.accept(TermFoldVisitor(self.log)) # get the atom
        return ast_goal # return the goal

    def visit_rule(self, ast_rule):
        ast_rule.head = ast_rule.head.accept(TermFoldVisitor(self.log)) # get the head
        ast_rule.consequence = ast_rule.consequence.accept(LogicalFoldVisitor(self.log)) # get the consequence
        return ast_rule # return the rule

    def visit_list(self, ast_list):
        term_visitor = TermFoldVisitor(self.log) # create a term visitor
        ast_list.terms = [term.accept(term_visitor) for term in ast_list.terms] # get the terms
        return ast_list # return the list

    def visit_literal(self, ast_literal):
        term_visitor = TermFoldVisitor(self.log) # create a term visitor
        ast_literal.terms = [term.accept(term_visitor) for term in ast_literal.terms] # get the terms
        ast_literal.annotations = [annotation.accept(term_visitor) for annotation in ast_literal.annotations] # get the annotations
        return ast_literal # return the literal


def validate(ast_agent, log):
    ast_agent = ast_agent.accept(ConstFoldVisitor(log)) # fold the constants

    for belief in ast_agent.beliefs: # iterate over the beliefs
        variables = list(belief.accept(FindVariablesVisitor())) # get the variables
        if variables: # if there are variables
            names = sorted(set(variable.name for variable in variables)) # get the names
            log.warning("implicit rule with unbound variables: %s (add ':- true' to acknowledge)",
                        ", ".join("'%s'" % name for name in names),
                        loc=belief.loc,
                        extra_locs=[variable.loc for variable in variables]) # log a warning

        for op in belief.accept(FindOpVisitor()): # iterate over the operators
            log.error("base belief can not contain this expression", loc=op.loc, extra_locs=[belief.loc])

    for rule in ast_agent.rules: # iterate over the rules
        for op in rule.head.accept(FindOpVisitor()): # iterate over the operators
            log.error("rule head is supposed to be unifiable, but contains non-const expression", loc=op.loc, extra_locs=[rule.loc])

    for plan in ast_agent.plans: # iterate over the plans
        for op in plan.event.head.accept(FindOpVisitor()): # iterate over the operators
            log.error("plan head is supposed to be unifiable, but contains non-const expression", loc=op.loc, extra_locs=[plan.loc])

        for annotation in plan.annotations: # iterate over the annotations
            log.warning("plan annotations are ignored as of yet", loc=annotation.loc, extra_locs=[plan.loc])

        if plan.event.goal_type != GoalType.belief and plan.event.trigger == Trigger.removal: # if the event is not a belief and the trigger is removal
            log.warning("recovery plans are ignored as of yet", loc=plan.loc)

    return ast_agent # return the agent


def parse(filename, tokens, log, included_files=frozenset(), directive=None):
    try: # try to
        return parse_agent(filename, tokens, log, included_files, directive=directive) # parse the agent
    except StopIteration: # if the iteration stops unexpectedly raise an error
        raise log.error("unexpected end of file", loc=tokens.peek() and tokens.peek().loc)


def main(source, hook):
    log = agentspeak.Log(agentspeak.get_logger(__name__), 3) # create a log object with a logger and a verbosity level

    tokens = agentspeak.lexer.TokenStream(source, log, 1) # create a token stream object with the source, the log and a verbosity level
    agent = parse(source.name, tokens, log) # parse the agent from the token stream object and the log object and get the agent

    log.throw() # throw the log object

    hook(agent) # call the hook with the agent 


def repl(hook):
    lineno = 0 # set the line number to 0
    tokens = [] # create a list of tokens

    while True:
        try:
            log = agentspeak.Log(agentspeak.get_logger(__name__), 3) # create a log object with a logger and a verbosity level

            if not tokens: # if there are no tokens
                line = agentspeak.util.prompt("agentspeak.parser >>> ") # get the line from the user with a prompt and a default value of "agentspeak.parser >>> "
            else: # if there are tokens
                line = agentspeak.util.prompt("agentspeak.parser ... ") # get the line from the user with a prompt and a default value of "agentspeak.parser ... "

            lineno += 1 # increment the line number

            tokens.extend(agentspeak.lexer.tokenize(agentspeak.StringSource("<stdin>", line), log, lineno)) # tokenize the line and add the tokens to the list of tokens

            while tokens: # while there are tokens in the list of tokens
                token_stream = iter(tokens) # create an iterator from the list of tokens and get the token stream object
                try:
                    agent = parse_agent("<stdin>", token_stream, log, frozenset()) # parse the agent from the token stream object and the log object and get the agent 
                except StopIteration:
                    log.throw() # throw the log object if there is an error
                    break
                else:
                    log.throw() 
                    hook(agent) 
                    tokens = list(token_stream) 
        except agentspeak.AggregatedError as error:
            print(str(error), file=sys.stderr) # print the error to stderr
            tokens = [] # reset the list of tokens
        except KeyboardInterrupt:
            print() # print a newline
            sys.exit(0) # exit with a status code of 0


if __name__ == "__main__":
    try:
        args = sys.argv[1:] # get the arguments from the command line and remove the first argument
        if args: # if there are arguments
            for arg in args: # iterate over the arguments
                with open(arg) as source: # open the file
                    main(source, print) # call the main function with the source and the print function
        elif sys.stdin.isatty(): # if there is no input from stdin
            repl(print) # call the repl function with the print function
        else:
            main(sys.stdin, print) # call the main function with the stdin and the print function
    except agentspeak.AggregatedError as error:
        print(str(error), file=sys.stderr) # print the error to stderr if there is an error 
        sys.exit(1) # exit with a status code of 1
